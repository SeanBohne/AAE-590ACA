function nu = kepler_true_anom(M, e)
% This function solves Kepler's equation to find true anomaly given mean
% anomaly

% Initial guess for the eccentric anomaly and iteration step size
E = M + e*sin(M)*(1 + e*cos(M));
delta_E = 1;
tol = 1e-8;

% Use Newton-Raphson method
while abs(delta_E) > tol
    f = E - e*sin(E) - M;
    fp = 1 - e*cos(E);
    delta_E = f/fp;
    E = E - delta_E;
end

% Calculate true anomaly from eccentric anomaly (rad)
nu = 2*atan2(sqrt(1 + e)*sin(E/2), sqrt(1 - e)*cos(E/2));
end