function dxdt = cartesian_pert(t, x, param)
% This function propogates cartesian states using two-body motion under
% various perturbations


dxdt = zeros(6, 1);
mu = param.mu_earth;

r = norm(x(1:3)); % Distance between orbited body and satellite (km)

% Store states to meaningful variables
x_pos = x(1);
y_pos = x(2);
z_pos = x(3);
r_vec = [x_pos;
         y_pos;
         z_pos];
x_vel = x(4);
y_vel = x(5);
z_vel = x(6);
v_vec = [x_vel;
         y_vel;
         z_vel];

dxdt(1:3) = [x_vel; 
             y_vel; 
             z_vel];
dxdt(4:6) = -mu*x(1:3)/(r^3);

if isfield(param, 'use_J2') && param.use_J2
    % Extract needed parameters
    r0 = param.r_earth;
    J2 = param.J2;
    [~, a_J2] = propagate.J2_pert(r_vec, v_vec, mu, r0, J2); % J2 perturbation accelerations (km/s^2)
    dxdt(4:6) = dxdt(4:6) + a_J2;
end

if isfield(param, 'use_SRP') && param.use_SRP
    % Extract needed parameters
    A_m = param.A_m;
    G0 = param.G0;
    omega_sun = param.omega_sun;
    AU = param.AU;
    [~, a_SRP] = propagate.SRP_pert(t, r_vec, v_vec, A_m, G0, omega_sun, AU); % SRP perturbation accelerations (km/s^2)
    dxdt(4:6) = dxdt(4:6) + a_SRP;
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
    dxdt(4:6) = dxdt(4:6) + a_Drag;
end
end