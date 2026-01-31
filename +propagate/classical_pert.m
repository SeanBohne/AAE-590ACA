function dxdt = classical_pert(t, x, param)
% This function propogates classical orbit elements under various 
% perturbations

dxdt = zeros(6, 1);
mu = param.mu_grav;

% Store states to meaningful variables
a = x(1); % Semi-major axis (km)
e = x(2); % Eccentricity
incl = x(3); % Inclination (rad)
RAAN = x(4); % Right ascension of ascending node (rad)
omega = x(5); % Argument of periapsis (rad)
M = x(6); % Mean anomaly (rad)
nu = elements.kepler_true_anom(M, e); % True anomaly (rad)

[r_vec, v_vec] = elements.classical_to_cartesian(a, e, incl, RAAN, omega, nu, mu); % Convert classical elements to ECI states

n = sqrt(mu/a^3); % Mean motion (rad/s)
p = a*(1 - e^2); % Semilatus rectum
b = a*sqrt(1 - e^2); % Perturbation parameter
r = p/(1 + e*cos(nu)); % Radial distance (km)
h = sqrt(mu*a*(1 - e^2)); % specific angular momentum(km^2/s)

B = 1/h*[2*a^2*e*sin(nu),             2*a^2*p/r,                0;
     p*sin(nu),                   (p + r)*cos(nu) + r*e,    0;
     0,                           0,                        r*cos(nu + omega);
     0,                           0,                        r*sin(nu + omega)/sin(incl);
     -p*cos(nu)/e,                (p + r)*sin(nu)/e,        -r*sin(nu + omega)/tan(incl);
     b*p*cos(nu)/(a*e) - 2*b*r/a, -b*(p + r)*sin(nu)/(a*e), 0]; % Perturbation mapping matrix

dxdt(6) = n;

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

if isfield(param, 'use_3body') && param.use_3body
    % Extract needed parameters
    GM_star = param.GM_star;
    l_star = param.l_star;
    mu_grav = param.mu_grav;
    [a_moon, ~] = propagate.moon_pert(t, r_vec, v_vec, GM_star, l_star, mu_grav); % Moon perturbation accelerations (km/s^2)
    dxdt(4:6) = dxdt(4:6) + a_moon;
end
end
