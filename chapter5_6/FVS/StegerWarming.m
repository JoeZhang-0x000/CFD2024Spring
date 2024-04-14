classdef StegerWarming < FVS
    % StegerWarming 分裂法
    methods(Static)
        function [f_p,f_n] = forward(U,gamma)
                        % 解析 U =[rho, rho * u, E]
            rho = U(1,:); % rho = rho
            u = U(2,:) ./ rho; % u = rho * u / rho
            p = (U(3,:) - 1/2 * rho .* (u).^2) * (gamma - 1);

            % 计算lambda_1, lambda_2, lambda_3
            num = length(u);
            c = sqrt(p .* gamma ./ rho);
            lam = zeros(3,num);
            lam(1,:) = u;
            lam(2,:) = u - c;
            lam(3,:) = u + c;
            lam_p = zeros(3,num); % postive lambda_k
            lam_n = zeros(3,num); % negative lambda_k
            
            lam_p = (lam + sqrt(lam.^2 + eps)) / 2;
            lam_n = (lam - sqrt(lam.^2 + eps)) / 2;

            [f_p, f_n] = FVS.computeFpFn(lam_p,lam_n,gamma,c,rho,u);
        end

    end
end