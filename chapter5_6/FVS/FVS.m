classdef FVS
    methods(Abstract)
        forward;
    end
    methods(Static)
        function [f_p,f_n] = computeFpFn(lam_p,lam_n,gamma,c,rho,u)
            % 计算f+, f-
            w = (3-gamma) * (lam_p(2,:)+lam_p(3,:)) .* (c.^2) / (gamma - 1) / 2;
            f_p = rho / (2*gamma) .* [
                2* (gamma-1) * lam_p(1,:) + lam_p(2,:) + lam_p(3,:);
                2* (gamma-1) * lam_p(1,:) .* u + lam_p(2,:) .* (u - c) + lam_p(3,:) .* (u + c);
                (gamma - 1) * lam_p(1,:) .* (u.^2) + lam_p(2,:) / 2 .* (u - c).^2 + lam_p(3,:) / 2 .* (u + c).^2 + w;
                ];

            w = (3-gamma) * (lam_n(2,:)+lam_n(3,:)) .* (c.^2) / (gamma - 1) / 2;
            f_n = rho / (2*gamma) .* [
                2* (gamma-1) * lam_n(1,:) + lam_n(2,:) + lam_n(3,:);
                2* (gamma-1) * lam_n(1,:) .* u + lam_n(2,:) .* (u - c) + lam_n(3,:) .* (u + c);
                (gamma - 1) * lam_n(1,:) .* (u).^2 + lam_n(2,:) / 2 .* (u - c).^2 + lam_n(3,:) / 2 .* (u + c).^2 + w;
                ];
        end
    end
end