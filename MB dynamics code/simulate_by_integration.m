function [t, y] = simulate_by_integration(func, tspan, y0)
    options = odeset('MaxStep', 0.05);

    [t, y] = ode113(func, tspan, y0, options);
end