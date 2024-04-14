classdef RK
    methods(Abstract)
       forward; % 用于管理Runge kutta方法的接口，所有的RK子类都实现该方法
    end
end