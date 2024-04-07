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


        function Q = stepFun(obj,FD,u)
            % 显式推进的方法
            % u_t = -f_x = Q 
            if ~ isa(FD,'FD')
                error('差分方法错误');
            end

            fvs = FVS;
            [f_p, f_n] = fvs.simple(u,obj.gamma);
            for i = 1:3
                tem = [f_p(i,:);f_n(i,:)];
                f_x(i,:) = FD.computeDerivative(tem,obj.dx);
            end
            
            Q = -f_x;
        end

        function obj = timeStep(obj,fd, rk)
            % 时间推进，默认三阶R-K
            func = @(u) obj.stepFun(fd,u);
            U = [obj.rho; obj.rho .* obj.u; obj.E]; % [3,num]
            U = rk.rk3(U,func,obj.dt); % 下一个时刻的 rho, rho * u, E

            % 更新
            obj.rho = U(1,:); % rho = rho
            obj.u = U(2,:) ./ obj.rho; % u = rho * u / rho
            obj.p = (U(3,:) - 1/2 * obj.rho .* (obj.u).^2) * (obj.gamma - 1);
            obj.t = obj.t + obj.dt; % 时间推进一步
            obj.update(); % 更新其他量
        end

    end
end

