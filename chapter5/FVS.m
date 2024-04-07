classdef FVS
    methods
        function [f_p,f_n] = sw(obj,U,gamma)
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
%             for i = 1:3
%                 for j = 1:num
%                     lam_p(i,j) = (lam(i,j) + sqrt(lam(i,j).^2+eps)) / 2;
%                     lam_n(i,j) = (lam(i,j) - sqrt(lam(i,j).^2+eps)) / 2;
%                 end
%             end
            lam_p = (lam + sqrt(lam.^2 + eps)) / 2;
            lam_n = (lam - sqrt(lam.^2 + eps)) / 2;


            % 计算f+, f-
            w = (3-gamma) * (lam_p(2,:)+lam_p(3,:)) * (c^2) / (gamma - 1) / 2;
            f_p = rho / (2*gamma) .* [
                2* (gamma-1) * lam_p(1,:) + lam_p(2,:) + lam_p(3,:);
                2* (gamma-1) * lam_p(1,:) .* u + lam_p(2,:) .* (u - c) + lam_p(3,:) .* (u + c);
                (gamma - 1) * lam_p(1,:) .* (u.^2) + lam_p(2,:) / 2 .* (u - c).^2 + lam_p(3,:) / 2 .* (u + c).^2 + w;
                ];

            w = (3-gamma) * (lam_n(2,:)+lam_n(3,:)) * (c^2) / (gamma - 1) / 2;
            f_n = rho / (2*gamma) .* [
                2* (gamma-1) * lam_n(1,:) + lam_n(2,:) + lam_n(3,:);
                2* (gamma-1) * lam_n(1,:) .* u + lam_n(2,:) .* (u - c) + lam_n(3,:) .* (u + c);
                (gamma - 1) * lam_n(1,:) .* (u).^2 + lam_n(2,:) / 2 .* (u - c).^2 + lam_n(3,:) / 2 .* (u + c).^2 + w;
                ];

        end

        function [f_p,f_n] = simple(obj,U,gamma)
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
            lam_p(lam>=0) = lam(lam>=0);
            lam_n(lam<0) = lam(lam<0);

            % 计算f+, f-
            w = (3-gamma) * (lam_p(2,:)+lam_p(3,:)) * (c^2) / (gamma - 1) / 2;
            f_p = rho / (2*gamma) .* [
                2* (gamma-1) * lam_p(1,:) + lam_p(2,:) + lam_p(3,:);
                2* (gamma-1) * lam_p(1,:) .* u + lam_p(2,:) .* (u - c) + lam_p(3,:) .* (u + c);
                (gamma - 1) * lam_p(1,:) .* (u.^2) + lam_p(2,:) / 2 .* (u - c).^2 + lam_p(3,:) / 2 .* (u + c).^2 + w;
                ];

            w = (3-gamma) * (lam_n(2,:)+lam_n(3,:)) * (c^2) / (gamma - 1) / 2;
            f_n = rho / (2*gamma) .* [
                2* (gamma-1) * lam_n(1,:) + lam_n(2,:) + lam_n(3,:);
                2* (gamma-1) * lam_n(1,:) .* u + lam_n(2,:) .* (u - c) + lam_n(3,:) .* (u + c);
                (gamma - 1) * lam_n(1,:) .* (u).^2 + lam_n(2,:) / 2 .* (u - c).^2 + lam_n(3,:) / 2 .* (u + c).^2 + w;
                ];
        end
    end
end