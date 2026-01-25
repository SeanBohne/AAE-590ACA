function dxdt = equinoctial_pert(t, x, param)
% This function propogates equinoctial orbit elements under various
% perturbations

dxdt = zeros(6, 1);
mu = param.mu_earth;

% Store states to meaningful variables
p = x(1); % Semi-major axis (km)
f = x(2); % Eccentricity
g = x(3); % Inclination (rad)
h = x(4); % Right ascension of ascending node (rad)
k = x(5); % Argument of periapsis (rad)
L = x(6); % Mean anomaly (rad)

[r_vec, v_vec] = elements.equinoctial_to_cartesian(p, f, g, h, k, L, mu); % Convert equinoctial elements to ECI states

q = 1 + f*cos(L) + g*sin(L); % ODE parameter
s_squared = 1 + h^2 + k^2; % ODE parameter

B = sqrt(p/mu)*[0,        2*p/q,                   0;
                sin(L),  ((q + 1)*cos(L) + f)/q, -g*(h*sin(L) - k*cos(L))/q;
                -cos(L), ((q + 1)*sin(L) + g)/q, f*(h*sin(L) - k*cos(L))/q;
                0,        0,                       s_squared*cos(L)/(2*q);
                0,        0,                       s_squared*sin(L)/(2*q);
                0,        0,                       (h*sin(L) - k*cos(L))/q]; % Perturbation mapping matrix

dxdt(6) = sqrt(mu*x(1))*(q/x(1))^2;

if isfield(param, 'use_J2') && param.use_J2
    % Extract needed parameters
    r0 = param.r_earth;
    J2 = param.J2;
    [a_J2, ~] = propagate.J2_pert(r_vec, v_vec, mu, r0, J2); % J2 perturbation accelerations (km/s^2)
    dxdt = dxdt + B*(a_J2);
end

if isfield(param, 'use_SRP') && param.use_SRP
    % Extract needed parameters
    A_m = param.A_m;
    G0 = param.G0;
    omega_sun = param.omega_sun;
    AU = param.AU;
    [a_SRP, ~] = propagate.SRP_pert(t, r_vec, v_vec, A_m, G0, omega_sun, AU); % SRP perturbation accelerations (km/s^2)
    dxdt = dxdt + B*(a_SRP);
end

if isfield(param, 'use_Drag') && param.use_Drag
    % Extract needed parameters
    A_m = param.A_m;
    rho_ref = param.rho_ref;
    dens_scale = param.dens_scale;
    h0 = param.h0;
    CD = param.CD;
    w_earth = param.w_earth;
    [a_Drag, ~] = propagate.Drag_pert(r_vec, v_vec, rho_ref, dens_scale, h0, CD, A_m, w_earth); % Drag perturbation accelerations (km/s^2)
    dxdt = dxdt + B*(a_Drag);
end
end