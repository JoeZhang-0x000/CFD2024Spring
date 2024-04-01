% FVS
classdef MySod < handle
    properties
        num; % 单元数量
        dx;
        x;
        x_min;
        x_max;
        u;
        c;
        rho;
        p;
        E;
        gamma;
        dt;
        t;
        ep = 1e-6;
        f_x;
    end

    methods
        function obj = MySod(dx,x_min,x_max,dt,gamma)
            obj.dx = dx;
            obj.x_min = x_min;
            obj.x_max = x_max;
            obj.gamma = gamma;
            obj.dt = dt;
            obj.t = 0;
        end

        function show(obj)
            subplot(311)
            plot(obj.x, obj.u)
            legend 'u'
            subplot(312)
            plot(obj.x, obj.rho)
            legend 'rho'
            subplot(313)
            plot(obj.x,obj.p)
            legend 'p'
        end

        function obj = build(obj,u_l,rho_l,p_l,u_r,rho_r,p_r)
            obj = mesh(obj);
            obj.u = ones(1,obj.num) * u_l;
            obj.u(obj.x > 0.5) = u_r;
            obj.rho = ones(1,obj.num) * rho_l;
            obj.rho(obj.x > 0.5) = rho_r;
            obj.p = ones(1,obj.num) * p_l;
            obj.p(obj.x > 0.5) = p_r;
            obj.update();
        end

        % 若修改了值，务必调用以更新其他相关的值
        function obj = update(obj)
            obj.c = sqrt(obj.gamma * obj.p / obj.rho);
            obj.E = obj.p / (obj.gamma - 1) + 1/2 * obj.rho .* (obj.u).^2;
        end

        function obj = mesh(obj)
            obj.x = obj.x_min : obj.dx : obj.x_max;
            obj.num = length(obj.x);
        end

        function fvs(obj)
            % 计算lambda_1, lambda_2, lambda_3
            lam = zeros(3,obj.num);
            lam(1,:) = obj.u;
            lam(2,:) = obj.u - obj.c;
            lam(3,:) = obj.u + obj.c;
            lam_p = zeros(3,obj.num); % postive lambda_k
            lam_n = zeros(3,obj.num); % negative lambda_k

            % 计算 lambda+, lambda-
            for i = 1:3
                lam_p(i,:) = (lam(i,:) + sqrt(lam(i,:).^2+obj.ep)) / 2;
                lam_n(i,:) = (lam(i,:) - sqrt(lam(i,:).^2+obj.ep)) / 2;
            end

            % 计算f+, f-
            w = (3-obj.gamma) * (lam_p(2,:)+lam_p(3,:)) .* (obj.c.^2) / (obj.gamma - 1) / 2;
            f_p = obj.rho / (2*obj.gamma) .* [
                2* (obj.gamma-1) * lam_p(1,:) + lam_p(2,:) + lam_p(3,:);
                2* (obj.gamma-1) * lam_p(1,:) .* obj.u + lam_p(2,:) .* (obj.u - obj.c) + lam_p(3,:) .* (obj.u + obj.c);
                (obj.gamma - 1) * lam_p(1,:) .* (obj.u).^2 + lam_p(2,:) / 2 .* (obj.u - obj.c).^2 + lam_p(3,:) / 2 .* (obj.u + obj.c).^2 + w;
                ];

            w = (3-obj.gamma) * (lam_n(2,:)+lam_n(3,:)) .* (obj.c.^2) / (obj.gamma - 1) / 2;
            f_n = obj.rho / (2*obj.gamma) .* [
                2* (obj.gamma-1) * lam_n(1,:) + lam_n(2,:) + lam_n(3,:);
                2* (obj.gamma-1) * lam_n(1,:) .* obj.u + lam_n(2,:) .* (obj.u - obj.c) + lam_n(3,:) .* (obj.u + obj.c);
                (obj.gamma - 1) * lam_n(1,:) .* (obj.u).^2 + lam_n(2,:) / 2 .* (obj.u - obj.c).^2 + lam_n(3,:) / 2 .* (obj.u + obj.c).^2 + w;
                ];
            %             disp(f_p);
            %             disp(f_n);
            % 利用不同迎风格式计算 \partial f^+ / \partial x, \partial f^- / \partial x
            % 默认调用nnd格式
            f_n_x = GVC2(f_n, obj.dx);
            f_p_x = GVC2(f_p, obj.dx);

            % 令 f_x = f_n_x + f_p_x
            obj.f_x = f_n_x + f_p_x;
        end

        function timeStep(obj, fvs, rk)
            % 时间推进，默认三阶R-K
            fd = @(u,dx) GVC2(u,dx);
            fun = @(u,dx) fvs.sw(u,dx,obj.gamma,obj.ep);
            U = [obj.rho; obj.rho .* obj.u; obj.E]; % [3,num]
            U = rk.rk3(U,fun,obj.dt,obj.dx); % 下一个时刻的 rho, rho * u, E

            % 更新
            obj.rho = U(1,:); % rho = rho
            obj.u = U(2,:) ./ obj.rho; % u = rho * u / rho
            obj.p = (U(3,:) - 1/2 * obj.rho .* (obj.u).^2) * (obj.gamma - 1);
            obj.t = obj.t + obj.dt; % 时间推进一步
            obj.update(); % 更新其他量
        end

    end
end

