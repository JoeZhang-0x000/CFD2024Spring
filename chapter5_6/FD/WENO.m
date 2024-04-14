classdef WENO < FD
    % 五阶WENO差分格式
    methods(Static)
        function Fx = forward(u,dx)
            [N,M] = size(u);
            if N ~= 2
                error('未定义');
            end
            F_p = u(1,:); % postive components
            F_n = u(2,:); % negative components
            Fh_p = zeros(1,M); % row 4 Fh_p_WENO, 正flux
            Fh_n = zeros(1,M); % 负flux
            Fh_p_c = zeros(3,M);
            Fh_n_c = zeros(3,M);
            Fx_p = zeros(1,M);
            Fx_n = zeros(1,M);
            Fx = zeros(1,M); % 差分结果
            IS_p = zeros(3,M);
            IS_n = zeros(3,M);
            Ome_p = zeros(3,M);
            Ome_n = zeros(3,M);
            Alpha_p = zeros(3,M);
            Alpha_n = zeros(3,M);
            p = 2;


            for j = 3:M-2
                Fh_p_c(1,j) = 1/3 * F_p(j-2) - 7/6 * F_p(j-1) + 11/6 * F_p(j);
                Fh_p_c(2,j) = -1/6 * F_p(j-1) + 5/6 * F_p(j) + 1/3 * F_p(j+1);
                Fh_p_c(3,j) = 1/3 * F_p(j) + 5/6 * F_p(j+1) - 1/6 * F_p(j+2);
                Fh_n_c(1,j) = 1/3 * F_n(j+2) - 7/6 * F_n(j+1) + 11/6 * F_n(j);
                Fh_n_c(2,j) = -1/6 * F_n(j+1) + 5/6 * F_n(j) + 1/3 * F_n(j-1);
                Fh_n_c(3,j) = 1/3 * F_n(j) + 5/6 * F_n(j-1) - 1/6 * F_n(j-2);

                IS_p(1,j) = 1/4 * (F_p(j-2) - 4 * F_p(j-1) + 3 * F_p(j))^2 + 13/12 * (F_p(j-2) - 2 * F_p(j-1) + F_p(j))^2;
                IS_p(2,j) = 1/4 * (F_p(j-1) - F_p(j+1))^2 + 13/12 * (F_p(j-1) -2 * F_p(j) + F_p(j+1))^2;
                IS_p(3,j) = 1/4 * (3* F_p(j) - 4 * F_p(j+1) + F_p(j+2))^2 + 13/12 * (F_p(j) - 2 * F_p(j+1) + F_p(j+2))^2;
                IS_n(1,j) = 1/4 * (F_n(j+2) - 4 * F_n(j+1) + 3 * F_n(j))^2 + 13/12 * (F_n(j+2) - 2 * F_n(j+1) + F_n(j))^2;
                IS_n(2,j) = 1/4 * (F_n(j+1) - F_n(j-1))^2 + 13/12 * (F_n(j+1) -2 * F_n(j) + F_n(j-1))^2;
                IS_n(3,j) = 1/4 * (3* F_n(j) - 4 * F_n(j-1) + F_n(j-2))^2 + 13/12 * (F_n(j) - 2 * F_n(j-1) + F_n(j-2))^2;

                Alpha_p(1,j) = 1/10 / (eps + IS_p(1,j))^p;
                Alpha_p(2,j) = 6/10 / (eps + IS_p(2,j))^p;
                Alpha_p(3,j) = 3/10 / (eps + IS_p(3,j))^p;
                Alpha_n(1,j) = 1/10 / (eps + IS_n(1,j))^p;
                Alpha_n(2,j) = 6/10 / (eps + IS_n(2,j))^p;
                Alpha_n(3,j) = 3/10 / (eps + IS_n(3,j))^p;

                Alpha_sum_p = sum(Alpha_p(:,j));
                Alpha_sum_n = sum(Alpha_n(:,j));

                for k = 1 : 3
                    Ome_p(k,j) = Alpha_p(k,j) / Alpha_sum_p;
                    Fh_p(j) = Fh_p(j) + Ome_p(k,j) * Fh_p_c(k,j);

                    Ome_n(k,j) = Alpha_n(k,j) / Alpha_sum_n;
                    Fh_n(j) = Fh_n(j) + Ome_n(k,j) * Fh_n_c(k,j);
                end
            end

            for j = 4 : M-3
                Fx_p(j) = (Fh_p(j) - Fh_p(j-1)) / dx;
                Fx_n(j) = (Fh_n(j+1) - Fh_n(j)) / dx;
                Fx(j) = Fx_p(j) + Fx_n(j);
            end

        end
    end
end