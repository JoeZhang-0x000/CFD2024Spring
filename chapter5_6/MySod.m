classdef MySod < handle
    properties
        num; % 单元数量
        dx;
        x;
        partition; % 隔板的位置(分界位置)
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
        discrete_method = WENO;
        time_forward_method = RK3;
        split_method = StegerWarming;
    end

    methods
        function obj = setDiscreteMethod(obj,new_method) % 设置空间离散方法
            if ~isa(new_method,'FD')
                error('离散方法不合法，仅支持继承自FD类')
            end
            obj.discrete_method = new_method;
        end
        
        function obj = setTimeForwardMethod(obj,new_method) % 设置时间推进方法
            if ~isa(new_method,'RK')
                error('时间推进方法不合法，仅支持继承自RK类')
            end
            obj.time_forward_method = new_method;
        end
        function obj = setSplitMethod(obj,new_method) % 设置分裂方法
            if ~isa(new_method,'FVS')
                error('分裂方法不合法，仅支持继承自FVS类')
            end
            obj.split_method = new_method;
        end

        function obj = MySod(dx,x_min,partition,x_max,dt,gamma)
            obj.dx = dx;
            obj.x_min = x_min;
            obj.partition = partition;
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

        function show3d(obj)
            colormap('jet')
            [xx,yy] = meshgrid(obj.x,[0,0.1]);
            subplot(311)
            contourf(xx,yy,repmat(obj.u,2,1));
            title 'u'
            colorbar;
            subplot(312)
            contourf(xx,yy,repmat(obj.rho,2,1));
            title 'rho'
            colorbar;
            subplot(313)
            contourf(xx,yy,repmat(obj.p,2,1));
            title 'p'
            colorbar;
        end

        function obj = build(obj,u_l,rho_l,p_l,u_r,rho_r,p_r)
            obj = mesh(obj);
            obj.u = ones(1,obj.num) * u_l;
            obj.u(obj.x >= obj.partition) = u_r;
            obj.rho = ones(1,obj.num) * rho_l;
            obj.rho(obj.x >= obj.partition) = rho_r;
            obj.p = ones(1,obj.num) * p_l;
            obj.p(obj.x >= obj.partition) = p_r;
            obj.update();
        end

        % 若修改了值，务必调用以更新其他相关的值
        function obj = update(obj)
            obj.c = sqrt(obj.gamma * obj.p / (obj.rho+eps));
            obj.E = obj.p / (obj.gamma - 1) + 1/2 * obj.rho .* (obj.u).^2;
        end

        function obj = mesh(obj)
            obj.x = obj.x_min : obj.dx : obj.x_max;
            obj.num = length(obj.x);
        end


        function Q = stepFun(obj,u)
            % 显式推进的方法
            % u_t = -f_x = Q 
            [f_p, f_n] = obj.split_method.forward(u,obj.gamma); % 调用FVS方法计算f^+,f^-
            f_x = zeros(3,obj.num); % 对f的每一项求导
            for i = 1:3
                tem = [f_p(i,:);f_n(i,:)];
                f_x(i,:) = obj.discrete_method.forward(tem,obj.dx);
            end
            
            Q = -f_x;
        end

        function obj = timeStep(obj)
            % 时间推进，默认三阶R-K
            func = @(u) obj.stepFun(u);
            U = [obj.rho; obj.rho .* obj.u; obj.E]; % [3,num]
            U = obj.time_forward_method.forward(U,func,obj.dt); % 下一个时刻的 rho, rho * u, E
    
            % 更新
            obj.rho = U(1,:); % rho = rho
            obj.u = U(2,:) ./ (obj.rho + eps); % u = rho * u / rho
            obj.p = (U(3,:) - 1/2 * obj.rho .* (obj.u).^2) * (obj.gamma - 1);
            obj.t = obj.t + obj.dt; % 时间推进一步
            obj.update(); % 更新其他量
        end

    end
end

