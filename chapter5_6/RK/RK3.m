classdef RK3 < RK
    methods(Static)
        function U = forward(U_0,Q,dt)
            U_1 = U_0 + dt * Q(U_0);
            U_2 = 3/4 * U_0 + 1/4 * U_1 + 1/4 * dt * Q(U_1);
            U = 1/3 * U_0 + 2/3 * U_2 + 2/3 * dt * Q(U_2);
        end
    end
end