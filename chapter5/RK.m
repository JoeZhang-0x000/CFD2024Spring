classdef RK
    methods
        function U = rk3(obj,U_0,Q,dt,dx)
            U_1 = U_0 + dt * Q(U_0,dx);
            U_2 = 3/4 * U_0 + 1/4 * U_1 + 1/4 * dt * Q(U_1,dx);
            U = 1/3 * U_0 + 2/3 * U_2 + 2/3 * dt * Q(U_2,dx);
        end
    end
end