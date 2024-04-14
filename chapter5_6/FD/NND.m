classdef NND < FD

    methods(Static)
        function Fx = forward(u,dx)
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
                Fh_p(j) = F_p(j) +     (0.5 * NND.Minmod((F_p(j) - F_p(j - 1)), (F_p(j + 1) - F_p(j))));
                Fh_n(j) = F_n(j + 1) - (0.5 * NND.Minmod((F_n(j + 1) - F_n(j)), (F_n(j + 2) - F_n(j + 1))));
                Fh(j) = Fh_p(j) + Fh_n(j);
            end

            xs_new = xs + 1;
            xt_new = xt;

            for j = xs_new : xt_new
                Fx(j) = (Fh(j) - Fh(j - 1)) / dx;
            end
        end
        function Q = Minmod(a, b)
            Q = 0;
            if a * b > 0
                Q = a;
                if abs(a) > abs(b)
                    Q = b;
                end
            end
        end
    end
end