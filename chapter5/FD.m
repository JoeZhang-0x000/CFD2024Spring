% % 差分方法的基类，创建时候指定不同的格式
classdef FD < handle
    properties
        id = 0;
    end

    methods
        % 根据id创建差分格式
        function obj = FD(id)
            obj.id = id;
        end

        function obj = setId(obj,id)
            obj.id = id;
        end

        % 根据method_id计算差分
        function u_x = computeDerivative(obj,u,dx)
            switch obj.id
                case 1 % 前差, (f_{j+1} - f_{j}) / dx
                    u_x = obj.forward(u,dx);
                case 2 % 后差, (f_{j} - f_{j-1}) / dx
                    u_x = obj.backward(u,dx);
                case 3 % TVD
                    u_x = obj.TVD(u,dx);
                case 4 % NNd
                    u_x = obj.NND(u,dx);
            end
        end

        function u_x = forward(obj,u,dx) % 前差, (f_{j+1} - f_{j}) / dx
            u_tem = [u(2:length(u)),u(length(u))]; % padding
            u_diff = u_tem - u;
            u_x = u_diff / dx;
        end

        function u_x = backward(obj,u,dx) % 后差, (f_{j} - f_{j-1}) / dx
            u_tem = [u(1),u(1:length(u)-1)]; % padding
            u_diff = u - u_tem;
            u_x = u_diff / dx;
        end

        function Fx = TVD(obj,u,dx) % TVD
            [N,M] = size(u);
            if N ~= 2
                error('未定义');
            end
            F_p = u(1,:);
            F_n = u(2,:);
            Fh_p = zeros(1,M);
            Fh_n = zeros(1,M);
            Fx_p = zeros(1,M);
            Fx_n = zeros(1,M);
            Fx = zeros(1,M);
            xs = 2;
            xt = M - xs;

            for j = xs : xt
                r_p = (F_p(j) - F_p(j - 1)) ./ (F_p(j + 1) - F_p(j) + eps); 
                r_n = (F_n(j + 2) - F_n(j + 1)) ./ (F_n(j + 1) - F_n(j) + eps);
                Phi_p = (r_p + abs(r_p)) ./ (1 + r_p);
                Phi_n = (r_n + abs(r_n)) ./ (1 + r_n);

                Fh_p(j) = F_p(j) + 0.5 * (Phi_p .* (F_p(j + 1) - F_p(j)));
                Fh_n(j) = F_n(j + 1) - 0.5 * (Phi_n .* (F_n(j + 1) - F_n(j)));
            end
            
            xs_new = xs + 1;
            xt_new = xt;

            for j = xs_new : xt_new
                Fx_p(j) = (Fh_p(j) - Fh_p(j - 1)) / dx;
                Fx_n(j) = (Fh_n(j) - Fh_n(j - 1)) / dx;
                Fx(j) = Fx_p(j) + Fx_n(j);
            end
        end

        function Q = Minmod(obj,a, b)
            Q = 0;
            if a * b > 0
                Q = a;
                if abs(a) > abs(b)
                    Q = b;
                end
            end

        end

        function Fx = NND(obj,u,dx) % NND
            [N,M] = size(u);
            if N ~= 2
                error('未定义');
            end
            F_p = u(1,:);
            F_n = u(2,:);
            Fh_p = zeros(1,M);
            Fh_n = zeros(1,M);
            Fx_p = zeros(1,M);
            Fx_n = zeros(1,M);
            Fh = zeros(1,M);
            Fx = zeros(1,M);
            xs = 2;
            xt = M - xs;

            for j = xs : xt
                Fh_p(j) = F_p(j) +     (0.5 * obj.Minmod((F_p(j) - F_p(j - 1)), (F_p(j + 1) - F_p(j))));
                Fh_n(j) = F_n(j + 1) - (0.5 * obj.Minmod((F_n(j + 1) - F_n(j)), (F_n(j + 2) - F_n(j + 1))));
                Fh(j) = Fh_p(j) + Fh_n(j);
            end

            xs_new = xs + 1;
            xt_new = xt;

            for j = xs_new : xt_new
                Fx(j) = (Fh(j) - Fh(j - 1)) / dx;
            end
        end

    end
end