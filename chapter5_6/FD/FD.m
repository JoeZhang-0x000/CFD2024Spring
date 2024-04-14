classdef FD < handle
    % 差分方法的基类，创建时候指定不同的格式
    methods(Abstract,Static)
        forward;
    end
end