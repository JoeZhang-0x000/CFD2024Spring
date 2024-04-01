classdef FVS
    methods
        function f_x = sw(obj,U,dx,gamma,ep)
            FD =@(u,dx) GVC2(u,dx);
            % 解析 U =[rho, rho * u, E]
            rho = U(1,:); % rho = rho
            u = U(2,:) ./ rho; % u = rho * u / rho
            p = (U(3,:) - 1/2 * rho .* (u).^2) * (gamma - 1);

            % 计算lambda_1, lambda_2, lambda_3
            num = length(u);
            c = sqrt(p * gamma / rho);
            lam = zeros(3,num);
            lam(1,:) = u;
            lam(2,:) = u - c;
            lam(3,:) = u + c;
            lam_p = zeros(3,num); % postive lambda_k
            lam_n = zeros(3,num); % negative lambda_k

            % 计算 lambda+, lambda-
            for i = 1:3
                lam_p(i,:) = (lam(i,:) + sqrt(lam(i,:).^2+ep)) / 2;
                lam_n(i,:) = (lam(i,:) - sqrt(lam(i,:).^2+ep)) / 2;
            end

            % 计算f+, f-
            w = (3-gamma) * (lam_p(2,:)+lam_p(3,:)) .* (c.^2) / (gamma - 1) / 2;
            f_p = rho / (2*gamma) .* [
                2* (gamma-1) * lam_p(1,:) + lam_p(2,:) + lam_p(3,:);
                2* (gamma-1) * lam_p(1,:) .* u + lam_p(2,:) .* (u - c) + lam_p(3,:) .* (u + c);
                (gamma - 1) * lam_p(1,:) .* (u).^2 + lam_p(2,:) / 2 .* (u - c).^2 + lam_p(3,:) / 2 .* (u + c).^2 + w;
                ];

            w = (3-gamma) * (lam_n(2,:)+lam_n(3,:)) .* (c.^2) / (gamma - 1) / 2;
            f_n = rho / (2*gamma) .* [
                2* (gamma-1) * lam_n(1,:) + lam_n(2,:) + lam_n(3,:);
                2* (gamma-1) * lam_n(1,:) .* u + lam_n(2,:) .* (u - c) + lam_n(3,:) .* (u + c);
                (gamma - 1) * lam_n(1,:) .* (u).^2 + lam_n(2,:) / 2 .* (u - c).^2 + lam_n(3,:) / 2 .* (u + c).^2 + w;
                ];

            % 利用不同迎风格式计算 \partial f^+ / \partial x, \partial f^- / \partial x
            f_n_x = FD(f_n, dx);
            f_p_x = FD(f_p, dx);
             
            f_x = f_n_x + f_p_x;

        end
    end
end