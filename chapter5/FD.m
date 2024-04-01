% % 差分方法的基类，创建时候指定不同的格式，调用时候使用forward()
% classdef FD < handle
%     properties
%         id = 0;
%     end
% 
%     methods 
%         % 根据id创建差分格式
%         function obj = FD(id)
%             obj.id = id;
%         end
% 
%         function obj = setId(obj,id)
%             obj.id = id;
%         end
%         
%         % 根据method_id计算差分
%         function u_x = forward(obj,u,dx)
%             n = length(u);
%             switch obj.id
%                 case -1 % 前差, (f_{j+1} - f_{j}) / dx
%                     u_t = [u(2:n)];
%                     u_diff = u_t - u(1:n-1);
%                     u_x = u_diff / dx;
%                 case 0 % 
%                     error('not implement error!');
%                 case 1 % 后差, (f_{j} - f_{j-1}) / dx
%                     u_t = [u(1:n-1)]; % 左移一位
%                     u_diff = u(2:n) - u_t;
%                     u_x = u_diff / dx;
%                 case 3 % GVC2
%                     % only cover 2 to n-1
%                     u_1 = [u(1:n-2)]; % u_{j-1} from 1 to n-2
%                     u_2 = [u(3:n)]; % u_{j+1}, from 3 to n
%                     u_ = [u(2:n-1)]; % u_{j}, from 2 to n-1
%                     limiter = (u_ - u_1) / (u_2 - u_);
%                     u_r = u_ + limiter .* (u_2 - u_) / 2; % u_{j+1/2}, from 2 to n-1
%                     u_l = u_1 + limiter .* 
%             end
%         end
% 
%     end
% end