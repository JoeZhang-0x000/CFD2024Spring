clear,clc;

templateFiled = struct('p',0,'u',0,'rho',0,'type','unkown');
numFiled = 5;
filed = repmat(templateFiled,numFiled,1);

filed(1).u = 0;
filed(1).rho = 1;
filed(1).p = 1;

filed(2).u = 0;
filed(2).rho = 0.125;
filed(2).p = 0.1;

gamma = 0.1;

% 第一步，判断左右两侧是激波还是膨胀波
[filed(1), filed(2)] = waveType(filed(1),filed(2),gamma);
fprintf('filed 1 %s\n',filed(1).type);
fprintf('filed 2 %s\n',filed(2).type);

% 第二步，求解中心区压力和速度
p_0 = (filed(1).p + filed(2).p)/2; % 初始预测的p*
f = @(p) F_1_2(p,filed(1),filed(2),gamma) - filed(1).u + filed(2).u;
p_0 = fzero(f,p_0); % 
u_0 = 1/2*(filed(1).u + filed(2).u + F(p_0,filed(2),gamma) - F(p_0,filed(1),gamma));
fprintf('中心区 p = %d, u = %d\n',p_0,u_0);

% 第三步，确定中心区接触间断两侧密度








% -------------------------FUNCTIONS---------------------------------------
% 判断给定区域的波形，rare or shock
function [filed_1, filed_2] = waveType(filed_1, filed_2,gamma)
    p_1 = filed_1.p;
    p_2 = filed_2.p;
    u_1 = filed_1.u;
    u_2 = filed_2.u;
    
    % 保证 p_1 >= p_2
    if p_2 > p_1
        disp('p2 > p1');
        [filed_1, filed_2] = waveType(filed_2, filed_1);
    end

    d_u = u_1 - u_2;

    if d_u >= F_1_2(p_1,filed_1,filed_2,gamma)
        % 情况一，左右都是激波
        filed_1.type = 'shock';
        filed_2.type = 'shock';
    elseif d_u >=F_1_2(p_2,filed_1,filed_2,gamma)
        % 情况二，左激波，右膨胀波
        filed_1.type = 'rarewave';
        filed_2.type = 'shock';
    elseif d_u >= F_1_2(0,filed_1,filed_2,gamma)
        % 情况四，左右都是膨胀波
        filed_1.type = 'rarewave';
        filed_2.type = 'rarewave';
    else
        % 情况五，真空
        filed_1.type = 'rarewave';
        filed_2.type = 'rarewave';
    end
end

% 根据给定的p*，以及1区或2区的状态，求f(p*,p_i,rho_i)
function f_p = F(p, filed, gamma)
    if p > filed.p
        f_p = F_1(p,filed,gamma);
    else
        f_p = F_2(p,filed,gamma);
    end
end

% f(p*,p_1,rho_1) + f(p*,p_2,rho_2)
function ret = F_1_2(p,filed_1,filed_2,gamma)
    ret = F(p,filed_1,gamma) + F(p,filed_2,gamma);
end

% p* > p_i
function ret = F_1(p_,filed, gamma)
    rho = filed.rho;
    p = filed.p;
    c = sqrt(gamma * p / rho);
    
    ret = (p_ - p)/(rho*c*((gamma+1)/(2*gamma)*(p_/p)+(gamma-1)/(2*gamma))^1/2);
end

% p* < p_i
function ret = F_2(p_,filed, gamma)
    rho = filed.rho;
    p = filed.p;
    c = sqrt(gamma * p / rho);

    ret = 2*c/(gamma-1)*((p_/p)^((gamma-1)/(2*gamma))-1);
end
