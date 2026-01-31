function dxdt = milankovitch_pert(t, x, param)
% This function propagates milankovitch orbit elements under various
% perturbations

dxdt = zeros(7, 1);
mu = param.mu_grav;

% Store states in meaningful variables
h_vec = x(1:3);
e_vec = x(4:6);
L = x(7);

h = norm(h_vec); % Angular momentum magnitude (km^2/s)
z_hat = [0, 0, 1]; % z direction unit vector

% Find radial distance (km)
[r_vec, v_vec] = elements.milankovitch_to_cartesian(h_vec, e_vec, L, mu);
r = norm(r_vec);

% Compute 3x3 cross product matrices
r_tilda = math.cross_mat(r_vec);
v_tilda = math.cross_mat(v_vec);
h_tilda = math.cross_mat(h_vec);

B = [r_tilda;
     1/mu*(v_tilda*r_tilda - h_tilda);
     dot(z_hat, r_vec)*h_vec'/(h*(h + dot(z_hat, h_vec)))]; % perturbation mapping matrix

dxdt(7) = h/r^2;

if isfield(param, 'use_J2') && param.use_J2
    % Extract needed parameters
    r0 = param.r_earth;
    J2 = param.J2;
    [~, a_J2] = propagate.J2_pert(r_vec, v_vec, mu, r0, J2); % J2 perturbation accelerations (km/s^2)
    dxdt = dxdt + B*(a_J2);
end

if isfield(param, 'use_SRP') && param.use_SRP
    % Extract needed parameters
    A_m = param.A_m;
    G0 = param.G0;
    omega_sun = param.omega_sun;
    AU = param.AU;
    [~, a_SRP] = propagate.SRP_pert(t, r_vec, v_vec, A_m, G0, omega_sun, AU); % SRP perturbation accelerations (km/s^2)
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
    [~, a_Drag] = propagate.Drag_pert(r_vec, v_vec, rho_ref, dens_scale, h0, CD, A_m, w_earth); % Drag perturbation accelerations (km/s^2)
    dxdt = dxdt + B*(a_Drag);
end

if isfield(param, 'use_3body') && param.use_3body
    % Extract needed parameters
    GM_star = param.GM_star;
    l_star = param.l_star;
    mu_grav = param.mu_grav;
    [~, a_moon] = propagate.moon_pert(t, r_vec, v_vec, GM_star, l_star, mu_grav); % Moon perturbation accelerations (km/s^2)
    dxdt(4:6) = dxdt(4:6) + a_moon;
end
end
