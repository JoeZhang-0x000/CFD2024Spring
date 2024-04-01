function u_x = GVC2(u,dx)
    [rn,cn] = size(u);
    u_x = zeros(rn,cn);
    % 扩充ghost cell
    U = [ones(rn,2).*u(:,1),u,ones(rn,1).*u(:,cn)];
    
    for i = 3:length(U)-1
        % 计算u_{j+1/2}
        u_r = U(:,i) + Limiter(U(:,i-1),U(:,i),U(:,i+1)) .* (U(:,i+1) - U(:,i)) / 2;
    
        % 计算u_{j-1/2}
        u_l = U(:,i-1) + Limiter(U(:,i-2),U(:,i-1),U(:,i)) .* (U(:,i) - U(:,i-1)) / 2;
    
        % 计算u_{j}_x
        u_x(:,i-2) = (u_r - u_l) / dx;

    end
end

% j-1, j, j+1
function ret = Limiter(u_1, u_2, u_3)
    [rn,cn] = size(u_1);
    ret = zeros(rn,1);
    t = (u_2  - u_1) ./ (u_3 - u_2);
    ret(t>=1) = 1;
    ret(t>0 & t<1) = t(t>0 & t<1);
end
