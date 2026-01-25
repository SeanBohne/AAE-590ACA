clc; clear; close all;

%% Part A
% Initialize simulation parameters
a = 9000; % Semi-major axis (km)
e = 0.2; % Eccentricity
incl = deg2rad(50); % Inclination (deg)
RAAN = deg2rad(10); % Right ascension of ascending node (deg)
omega = deg2rad(20); % Argument of periapsis (deg)
nu = deg2rad(0); % True anomaly (deg)
M = elements.kepler_mean_anom(nu, e); % Mean anomaly (deg)

param.mu_earth = 3.986*10^5; % Earth gravitational parameter (km^3/s^2)
param.r_earth = 6378.1; % Earth radius (km)
param.w_earth = 7.292115*10^-5; % Earth rotation rate (rad/s)
param.omega_sun = 2*pi/(365.25*24*3600); % Earth rotation around sun (rad/s)
param.A_m = 5.4e-6; % Area-to-mass ratio (km^2/kg)
param.G0 = 1.02e14; % Solar flux constant (kg*km/s^2)
param.AU = 149597870.7; % Astronomical unit (km)
param.J2 = 1.0826e-3; % J2 perturbation parameter
param.rho_ref = 1.225e9; % Reference density (kg/km^3)
param.dens_scale = 1/7.8; % Density scaling parameter (1/km)
param.h0 = param.r_earth; % reference altitude (km)
param.CD = 2.1; % Drag coefficient

num_orbits = 15; % number of simulated orbits
T = num_orbits*2*pi*sqrt(a^3/param.mu_earth); % Orbital period (s)
t = linspace(0, T, 10*T); % Time vector (s)

% Test Kepler's equation solver
nu_test = deg2rad(20); % Test true anomaly (deg)
M_test = elements.kepler_mean_anom(nu_test, e);
disp('Mean anomaly test value: ')
disp(rad2deg(M_test))

nu_test1 = elements.kepler_true_anom(M_test, e);
disp('True anomaly test value: ')
disp(rad2deg(nu_test1))

% Test element conversions to ECI
[r_eci, v_eci] = elements.classical_to_cartesian(a, e, incl, RAAN, omega, nu, param.mu_earth);
[a1, e1, incl1, RAAN1, omega1, nu1] = elements.cartesian_to_classical(r_eci, v_eci, param.mu_earth);
classical_residual = [a; e; incl; RAAN; omega; nu] - [a1; e1; incl1; RAAN1; omega1; nu1];

disp('Classical residuals: ')
disp(classical_residual)

[p, f, g, h, k, L] = elements.classical_to_equinoctial(a, e, incl, RAAN, omega, nu);
[r_eci1, v_eci1] = elements.equinoctial_to_cartesian(p, f, g, h, k, L, param.mu_earth);
[a2, e2, incl2, RAAN2, omega2, nu2] = elements.cartesian_to_classical(r_eci1, v_eci1, param.mu_earth);
equinoctial_residual = [a; e; incl; RAAN; omega; nu] - [a2; e2; incl2; RAAN2; omega2; nu2];

disp('Equinoctial residuals: ')
disp(equinoctial_residual)

[h_vec, e_vec, L] = elements.classical_to_milankovitch(a, e, incl, RAAN, omega, nu, param.mu_earth);
[r_eci2, v_eci2] = elements.milankovitch_to_cartesian(h_vec, e_vec, L, param.mu_earth);
[a3, e3, incl3, RAAN3, omega3, nu3] = elements.cartesian_to_classical(r_eci2, v_eci2, param.mu_earth);
milankovitch_residual = [a; e; incl; RAAN; omega; nu] - [a3; e3; incl3; RAAN3; omega3; nu3];

disp('Milankovitch residuals: ')
disp(milankovitch_residual)

%% Part B
% Define variables to specify perturbations used in simulation. false =
% excluded, true = included
param.use_J2 = false;
param.use_SRP = false;
param.use_Drag = false;

% Simulation parameters
opts1 = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); % ODE solver options 1
opts2 = odeset('RelTol', 1e-9, 'AbsTol', 1e-9); % ODE solver options 2
opts3 = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); % ODE solver options 3
opts4 = odeset('RelTol', 1e-3, 'AbsTol', 1e-3); % ODE solver options 4
opts = {opts1, opts2, opts3, opts4}; % Combine ODE solver options

% Preallocate variables
Nt = length(t);
sol_cart_b = zeros(Nt, 6, length(opts));
sol_eci_classical_b = zeros(6, Nt, length(opts));
sol_eci_equinoctial_b = zeros(6, Nt, length(opts));
sol_eci_milankovitch_b = zeros(6, Nt, length(opts));

for z = 1:length(opts)
    opts_z = opts{z};

    % Simulate trajectory using cartesian state
    x0_cart = [r_eci; v_eci];
    [t_cart_b, temp_cart_b] = ode45(@propagate.cartesian, t, x0_cart, opts_z, param.mu_earth);
    sol_cart_b(:, :, z) = temp_cart_b;

    % Simulate trajectory using keplerian orbit elements
    x0_classical_b = [a; e; incl; RAAN; omega; M];
    [t_classical_b, sol_classical_b] = ode45(@propagate.classical, t, x0_classical_b, opts_z, param.mu_earth);
    sol_classical_b = sol_classical_b';

    % Preallocate ECI position and velocity vectors
    r_eci_classical_b = zeros(length(x0_classical_b)/2, length(sol_classical_b));
    v_eci_classical_b = zeros(length(x0_classical_b)/2, length(sol_classical_b));

    % Convert classical orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_classical_b)
        a_sol = sol_classical_b(1, i);
        e_sol = sol_classical_b(2, i);
        incl_sol = sol_classical_b(3, i);
        RAAN_sol = sol_classical_b(4, i);
        omega_sol = sol_classical_b(5, i);
        nu_sol = elements.kepler_true_anom(sol_classical_b(6, i), e_sol);
        [r_eci_classical_b(:, i), v_eci_classical_b(:, i)] = elements.classical_to_cartesian(a_sol, e_sol, incl_sol, RAAN_sol, omega_sol, nu_sol, param.mu_earth);
    end
    sol_eci_classical_b(:, :, z) = [r_eci_classical_b; v_eci_classical_b]; % Combine cartesian state vector in ECI

    % Simulate trajectory using equinoctial orbit elements
    [p, f, g, h, k, L] = elements.classical_to_equinoctial(a, e, incl, RAAN, omega, nu);
    x0_equinoctial = [p, f, g, h, k, L];
    [t_equinoctial_b, sol_equinoctial_b] = ode45(@propagate.equinoctial, t, x0_equinoctial, opts_z, param.mu_earth);
    sol_equinoctial_b = sol_equinoctial_b';

    % Preallocate ECI position and velocity vectors
    r_eci_equinoctial_b = zeros(length(x0_equinoctial)/2, length(sol_equinoctial_b));
    v_eci_equinoctial_b = zeros(length(x0_equinoctial)/2, length(sol_equinoctial_b));

    % Convert equinoctial orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_equinoctial_b)
        p_sol = sol_equinoctial_b(1, i);
        f_sol = sol_equinoctial_b(2, i);
        g_sol = sol_equinoctial_b(3, i);
        h_sol = sol_equinoctial_b(4, i);
        k_sol = sol_equinoctial_b(5, i);
        L_sol = sol_equinoctial_b(6, i);
        [r_eci_equinoctial_b(:, i), v_eci_equinoctial_b(:, i)] = elements.equinoctial_to_cartesian(p_sol, f_sol, g_sol, h_sol, k_sol, L_sol, param.mu_earth);
    end
    sol_eci_equinoctial_b(:, :, z) = [r_eci_equinoctial_b; v_eci_equinoctial_b]; % Combine cartesian state vector in ECI

    % Simulate trajectory using milankovitch orbit elements
    [h_vec, e_vec, L] = elements.classical_to_milankovitch(a, e, incl, RAAN, omega, nu, param.mu_earth);
    x0_milankovitch_b = [h_vec; e_vec; L];
    [t_milankovitch_b, sol_milankovitch_b] = ode45(@propagate.milankovitch, t, x0_milankovitch_b, opts_z, param.mu_earth);
    sol_milankovitch_b = sol_milankovitch_b';

    h_sol_milankovitch_b = sol_milankovitch_b(1:3, :); % Simulate milankovitch angular momentum vector (km^2/s)
    e_sol_milankovitch_b = sol_milankovitch_b(4:6, :); % Simulate milankovitch eccentricity vector

    milankovitch_constraint_b = zeros(1, length(sol_milankovitch_b), length(opts)); % Preallocate constraint variable

    % Calculate milankovitch constraint over time
    for i = 1:length(sol_milankovitch_b)
        milankovitch_constraint_b(1, i, z) = dot(h_sol_milankovitch_b(:, i), e_sol_milankovitch_b(:, i));
    end

    % Preallocate ECI position and velocity vectors
    r_eci_milankovitch_b = zeros(floor(length(x0_milankovitch_b)/2), length(sol_milankovitch_b));
    v_eci_milankovitch_b = zeros(floor(length(x0_milankovitch_b)/2), length(sol_milankovitch_b));

    % Convert equinoctial orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_milankovitch_b)
        h_sol_vec = sol_milankovitch_b(1:3, i);
        e_sol_vec = sol_milankovitch_b(4:6, i);
        L_sol = sol_milankovitch_b(7, i);
        [r_eci_milankovitch_b(:, i), v_eci_milankovitch_b(:, i)] = elements.milankovitch_to_cartesian(h_sol_vec, e_sol_vec, L_sol, param.mu_earth);
    end
    sol_eci_milankovitch_b(:, :, z) = [r_eci_milankovitch_b; v_eci_milankovitch_b]; % Combine cartesian state vector in ECI
end

% Plot orbit simulations on one plot
for z = 1:length(opts)
    figure(z)
    plot3(sol_cart_b(:, 1, z), sol_cart_b(:, 2, z), sol_cart_b(:, 3, z))
    hold on
    plot3(sol_eci_classical_b(1, :, z), sol_eci_classical_b(2, :, z), sol_eci_classical_b(3, :, z))
    plot3(sol_eci_equinoctial_b(1, :, z), sol_eci_equinoctial_b(2, :, z), sol_eci_equinoctial_b(3, :, z))
    plot3(sol_eci_milankovitch_b(1, :, z), sol_eci_milankovitch_b(2, :, z), sol_eci_milankovitch_b(3, :, z))
    title('Unperturbed Orbit Trajectory - Sean Bohne', 'FontSize', 8)
    xlabel('x Position (km)')
    ylabel('y Position (km)')
    zlabel('z Position (km)')
    legend('Cartesian', 'Keplerian', 'Equinoctial', 'Milankovitch', 'location', 'best')
    axis equal
    grid on
end

% Plot milankovitch constraint over time
figure
for z = 1:length(opts)
    plot(t_milankovitch_b/3600, milankovitch_constraint_b(1, :, z))
    hold on
end
title('Milankovitch Constraint Over Time - Sean Bohne', 'FontSize', 8)
xlabel('Time (hrs)')
ylabel('Milankovitch Constraint')
grid on
legend('1e-12', '1e-9', '1e-6', '1e-3', 'location', 'best')

%% Part C
% Define variables to specify perturbations used in simulation. false =
% excluded, true = included
param.use_J2 = true;
param.use_SRP = false;
param.use_Drag = false;

% Preallocate variables
Nt = length(t);
sol_cart_b = zeros(Nt, 6, length(opts));
sol_eci_classical_b = zeros(6, Nt, length(opts));
sol_eci_equinoctial_b = zeros(6, Nt, length(opts));
sol_eci_milankovitch_b = zeros(6, Nt, length(opts));

for z = 1:length(opts)
    opts_z = opts{z};
    % Simulate trajectory under J2 perturbations using cartesian state
    x0_cart = [r_eci; v_eci];
    [t_cart_c, temp_cart_c] = ode45(@propagate.cartesian_pert, t, x0_cart, opts_z, param);
    sol_cart_c(:, :, z) = temp_cart_c;

    % Simulate trajectory using keplerian orbit elements
    x0_classical_c = [a; e; incl; RAAN; omega; M];
    [t_classical_c, sol_classical_c] = ode45(@propagate.classical_pert, t, x0_classical_c, opts_z, param);
    sol_classical_c = sol_classical_c';

    % Preallocate ECI position and velocity vectors
    r_eci_classical_c = zeros(length(x0_classical_c)/2, length(sol_classical_c));
    v_eci_classical_c = zeros(length(x0_classical_c)/2, length(sol_classical_c));

    % Convert classical orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_classical_c)
        a_sol = sol_classical_c(1, i);
        e_sol = sol_classical_c(2, i);
        incl_sol = sol_classical_c(3, i);
        RAAN_sol = sol_classical_c(4, i);
        omega_sol = sol_classical_c(5, i);
        nu_sol = elements.kepler_true_anom(sol_classical_c(6, i), e_sol);
        [r_eci_classical_c(:, i), v_eci_classical_c(:, i)] = elements.classical_to_cartesian(a_sol, e_sol, incl_sol, RAAN_sol, omega_sol, nu_sol, param.mu_earth);
    end
    sol_eci_classical_c(:, :, z) = [r_eci_classical_c; v_eci_classical_c]; % Combine cartesian state vector in ECI

    % Simulate trajectory using equinoctial orbit elements
    [p, f, g, h, k, L] = elements.classical_to_equinoctial(a, e, incl, RAAN, omega, nu);
    x0_equinoctial = [p, f, g, h, k, L];
    [t_equinoctial_c, sol_equinoctial_c] = ode45(@propagate.equinoctial_pert, t, x0_equinoctial, opts_z, param);
    sol_equinoctial_c = sol_equinoctial_c';

    % Preallocate ECI position and velocity vectors
    r_eci_equinoctial_c = zeros(length(x0_equinoctial)/2, length(sol_equinoctial_c));
    v_eci_equinoctial_c = zeros(length(x0_equinoctial)/2, length(sol_equinoctial_c));

    % Convert equinoctial orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_equinoctial_c)
        p_sol = sol_equinoctial_c(1, i);
        f_sol = sol_equinoctial_c(2, i);
        g_sol = sol_equinoctial_c(3, i);
        h_sol = sol_equinoctial_c(4, i);
        k_sol = sol_equinoctial_c(5, i);
        L_sol = sol_equinoctial_c(6, i);
        [r_eci_equinoctial_c(:, i), v_eci_equinoctial_c(:, i)] = elements.equinoctial_to_cartesian(p_sol, f_sol, g_sol, h_sol, k_sol, L_sol, param.mu_earth);
    end
    sol_eci_equinoctial_c(:, :, z) = [r_eci_equinoctial_c; v_eci_equinoctial_c]; % Combine cartesian state vector in ECI

    % Simulate trajectory using milankovitch orbit elements
    [h_vec, e_vec, L] = elements.classical_to_milankovitch(a, e, incl, RAAN, omega, nu, param.mu_earth);
    x0_milankovitch_c = [h_vec; e_vec; L];
    [t_milankovitch_c, sol_milankovitch_c] = ode45(@propagate.milankovitch_pert, t, x0_milankovitch_c, opts_z, param);
    sol_milankovitch_c = sol_milankovitch_c';

    h_sol_milankovitch_c = sol_milankovitch_c(1:3, :); % Simulate milankovitch angular momentum vector (km^2/s)
    e_sol_milankovitch_c = sol_milankovitch_c(4:6, :); % Simulate milankovitch eccentricity vector

    milankovitch_constraint_c = zeros(1, length(sol_milankovitch_c), length(opts)); % Preallocate constraint variable

    % Calculate milankovitch constraint over time
    for i = 1:length(sol_milankovitch_c)
        milankovitch_constraint_c(1, i, z) = dot(h_sol_milankovitch_c(:, i), e_sol_milankovitch_c(:, i));
    end

    % Preallocate ECI position and velocity vectors
    r_eci_milankovitch_c = zeros(floor(length(x0_milankovitch_c)/2), length(sol_milankovitch_c));
    v_eci_milankovitch_c = zeros(floor(length(x0_milankovitch_c)/2), length(sol_milankovitch_c));

    % Convert equinoctial orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_milankovitch_c)
        h_sol_vec = sol_milankovitch_c(1:3, i);
        e_sol_vec = sol_milankovitch_c(4:6, i);
        L_sol = sol_milankovitch_c(7, i);
        [r_eci_milankovitch_c(:, i), v_eci_milankovitch_c(:, i)] = elements.milankovitch_to_cartesian(h_sol_vec, e_sol_vec, L_sol, param.mu_earth);
    end
    sol_eci_milankovitch_c(:, :, z) = [r_eci_milankovitch_c; v_eci_milankovitch_c]; % Combine cartesian state vector in ECI
end

% Plot orbit simulations on one plot

for z = 1:length(opts)
    figure(z + length(opts) + 1)
    plot3(sol_cart_c(:, 1, z), sol_cart_c(:, 2, z), sol_cart_c(:, 3, z))
    hold on
    plot3(sol_eci_classical_c(1, :, z), sol_eci_classical_c(2, :, z), sol_eci_classical_c(3, :, z))
    plot3(sol_eci_equinoctial_c(1, :, z), sol_eci_equinoctial_c(2, :, z), sol_eci_equinoctial_c(3, :, z))
    plot3(sol_eci_milankovitch_c(1, :, z), sol_eci_milankovitch_c(2, :, z), sol_eci_milankovitch_c(3, :, z))
    title('Orbit Simulation Under Earth J2 Perturbations - Sean Bohne', 'FontSize', 8)
    xlabel('x Position (km)')
    ylabel('y Position (km)')
    zlabel('z Position (km)')
    legend('Cartesian', 'Keplerian', 'Equinoctial', 'Milankovitch', 'location', 'best')
    axis equal
    grid on
end

% Plot milankovitch constraint over time
figure
for z = 1:length(opts)
    plot(t_milankovitch_c/3600, milankovitch_constraint_c(1, :, z))
    hold on
end
title('Milankovitch Constraint Over Time - Sean Bohne', 'FontSize', 8)
xlabel('Time (hrs)')
ylabel('Milankovitch Constraint')
grid on
legend('1e-12', '1e-9', '1e-6', '1e-3', 'location', 'southeast')

%% Part D
% Define variables to specify perturbations used in simulation. false =
% excluded, true = included
param.use_J2 = true;
param.use_SRP = true;
param.use_Drag = false;
% Preallocate variables
Nt = length(t);
sol_cart_d = zeros(Nt, 6, length(opts));
sol_eci_classical_d = zeros(6, Nt, length(opts));
sol_eci_equinoctial_d = zeros(6, Nt, length(opts));
sol_eci_milankovitch_d = zeros(6, Nt, length(opts));

for z = 1:length(opts)
    opts_z = opts{z};
    % Simulate trajectory under J2 perturbations using cartesian state
    x0_cart = [r_eci; v_eci];
    [t_cart_d, temp_cart_d] = ode45(@propagate.cartesian_pert, t, x0_cart, opts_z, param);
    sol_cart_d(:, :, z) = temp_cart_d;

    % Simulate trajectory using keplerian orbit elements
    x0_classical_c = [a; e; incl; RAAN; omega; M];
    [t_classical_d, sol_classical_d] = ode45(@propagate.classical_pert, t, x0_classical_c, opts_z, param);
    sol_classical_d = sol_classical_d';

    % Preallocate ECI position and velocity vectors
    r_eci_classical_d = zeros(length(x0_classical_c)/2, length(sol_classical_d));
    v_eci_classical_d = zeros(length(x0_classical_c)/2, length(sol_classical_d));

    % Convert classical orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_classical_d)
        a_sol = sol_classical_d(1, i);
        e_sol = sol_classical_d(2, i);
        incl_sol = sol_classical_d(3, i);
        RAAN_sol = sol_classical_d(4, i);
        omega_sol = sol_classical_d(5, i);
        nu_sol = elements.kepler_true_anom(sol_classical_d(6, i), e_sol);
        [r_eci_classical_d(:, i), v_eci_classical_d(:, i)] = elements.classical_to_cartesian(a_sol, e_sol, incl_sol, RAAN_sol, omega_sol, nu_sol, param.mu_earth);
    end
    sol_eci_classical_d(:, :, z) = [r_eci_classical_d; v_eci_classical_d]; % Combine cartesian state vector in ECI

    % Simulate trajectory using equinoctial orbit elements
    [p, f, g, h, k, L] = elements.classical_to_equinoctial(a, e, incl, RAAN, omega, nu);
    x0_equinoctial = [p, f, g, h, k, L];
    [t_equinoctial_d, sol_equinoctial_d] = ode45(@propagate.equinoctial_pert, t, x0_equinoctial, opts_z, param);
    sol_equinoctial_d = sol_equinoctial_d';

    % Preallocate ECI position and velocity vectors
    r_eci_equinoctial_d = zeros(length(x0_equinoctial)/2, length(sol_equinoctial_d));
    v_eci_equinoctial_d = zeros(length(x0_equinoctial)/2, length(sol_equinoctial_d));

    % Convert equinoctial orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_equinoctial_d)
        p_sol = sol_equinoctial_d(1, i);
        f_sol = sol_equinoctial_d(2, i);
        g_sol = sol_equinoctial_d(3, i);
        h_sol = sol_equinoctial_d(4, i);
        k_sol = sol_equinoctial_d(5, i);
        L_sol = sol_equinoctial_d(6, i);
        [r_eci_equinoctial_d(:, i), v_eci_equinoctial_d(:, i)] = elements.equinoctial_to_cartesian(p_sol, f_sol, g_sol, h_sol, k_sol, L_sol, param.mu_earth);
    end
    sol_eci_equinoctial_d(:, :, z) = [r_eci_equinoctial_d; v_eci_equinoctial_d]; % Combine cartesian state vector in ECI

    % Simulate trajectory using milankovitch orbit elements
    [h_vec, e_vec, L] = elements.classical_to_milankovitch(a, e, incl, RAAN, omega, nu, param.mu_earth);
    x0_milankovitch_d = [h_vec; e_vec; L];
    [t_milankovitch_d, sol_milankovitch_d] = ode45(@propagate.milankovitch_pert, t, x0_milankovitch_d, opts_z, param);
    sol_milankovitch_d = sol_milankovitch_d';

    h_sol_milankovitch_d = sol_milankovitch_d(1:3, :); % Simulate milankovitch angular momentum vector (km^2/s)
    e_sol_milankovitch_d = sol_milankovitch_d(4:6, :); % Simulate milankovitch eccentricity vector

    milankovitch_constraint_d = zeros(1, length(sol_milankovitch_d), length(opts)); % Preallocate constraint variable

    % Calculate milankovitch constraint over time
    for i = 1:length(sol_milankovitch_d)
        milankovitch_constraint_d(1, i, z) = dot(h_sol_milankovitch_d(:, i), e_sol_milankovitch_d(:, i));
    end

    % Preallocate ECI position and velocity vectors
    r_eci_milankovitch_d = zeros(floor(length(x0_milankovitch_d)/2), length(sol_milankovitch_d));
    v_eci_milankovitch_d = zeros(floor(length(x0_milankovitch_d)/2), length(sol_milankovitch_d));

    % Convert equinoctial orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_milankovitch_d)
        h_sol_vec = sol_milankovitch_d(1:3, i);
        e_sol_vec = sol_milankovitch_d(4:6, i);
        L_sol = sol_milankovitch_d(7, i);
        [r_eci_milankovitch_d(:, i), v_eci_milankovitch_d(:, i)] = elements.milankovitch_to_cartesian(h_sol_vec, e_sol_vec, L_sol, param.mu_earth);
    end
    sol_eci_milankovitch_d(:, :, z) = [r_eci_milankovitch_d; v_eci_milankovitch_d]; % Combine cartesian state vector in ECI
end

% Plot orbit simulations on one plot

for z = 1:length(opts)
    figure(z + 2*length(opts) + 2)
    plot3(sol_cart_d(:, 1, z), sol_cart_d(:, 2, z), sol_cart_d(:, 3, z))
    hold on
    plot3(sol_eci_classical_d(1, :, z), sol_eci_classical_d(2, :, z), sol_eci_classical_d(3, :, z))
    plot3(sol_eci_equinoctial_d(1, :, z), sol_eci_equinoctial_d(2, :, z), sol_eci_equinoctial_d(3, :, z))
    plot3(sol_eci_milankovitch_d(1, :, z), sol_eci_milankovitch_d(2, :, z), sol_eci_milankovitch_d(3, :, z))
    title('Orbit Simulation Under Earth J2 and SRP Perturbations - Sean Bohne', 'FontSize', 8)
    xlabel('x Position (km)')
    ylabel('y Position (km)')
    zlabel('z Position (km)')
    legend('Cartesian', 'Keplerian', 'Equinoctial', 'Milankovitch', 'location', 'best')
    axis equal
    grid on
end

% Plot milankovitch constraint over time
figure
for z = 1:length(opts)
    plot(t_milankovitch_d/3600, milankovitch_constraint_d(1, :, z))
    hold on
end
title('Milankovitch Constraint Over Time - Sean Bohne', 'FontSize', 8)
xlabel('Time (hrs)')
ylabel('Milankovitch Constraint')
grid on
legend('1e-12', '1e-9', '1e-6', '1e-3', 'location', 'southeast')

%% Part E
% Define variables to specify perturbations used in simulation. false =
% excluded, true = included
param.use_J2 = true;
param.use_SRP = true;
param.use_Drag = true;

% Preallocate variables
Nt = length(t);
sol_cart_d = zeros(Nt, 6, length(opts));
sol_eci_classical_e = zeros(6, Nt, length(opts));
sol_eci_equinoctial_e = zeros(6, Nt, length(opts));
sol_eci_milankovitch_e = zeros(6, Nt, length(opts));

for z = 1:length(opts)
    opts_z = opts{z};
    % Simulate trajectory under J2 perturbations using cartesian state
    x0_cart = [r_eci; v_eci];
    [t_cart_e, temp_cart_e] = ode45(@propagate.cartesian_pert, t, x0_cart, opts_z, param);
    sol_cart_d(:, :, z) = temp_cart_e;

    % Simulate trajectory using keplerian orbit elements
    x0_classical_e = [a; e; incl; RAAN; omega; M];
    [t_classical_e, sol_classical_e] = ode45(@propagate.classical_pert, t, x0_classical_e, opts_z, param);
    sol_classical_e = sol_classical_e';
    
    % Preallocate ECI position and velocity vectors
    r_eci_classical_e = zeros(length(x0_classical_e)/2, length(sol_classical_e));
    v_eci_classical_e = zeros(length(x0_classical_e)/2, length(sol_classical_e));
    
    % Convert classical orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_classical_e)
        a_sol = sol_classical_e(1, i);
        e_sol = sol_classical_e(2, i);
        incl_sol = sol_classical_e(3, i);
        RAAN_sol = sol_classical_e(4, i);
        omega_sol = sol_classical_e(5, i);
        nu_sol = elements.kepler_true_anom(sol_classical_e(6, i), e_sol);
        [r_eci_classical_e(:, i), v_eci_classical_e(:, i)] = elements.classical_to_cartesian(a_sol, e_sol, incl_sol, RAAN_sol, omega_sol, nu_sol, param.mu_earth);
    end
    sol_eci_classical_e(:, :, z) = [r_eci_classical_e; v_eci_classical_e]; % Combine cartesian state vector in ECI
    
    % Simulate trajectory using equinoctial orbit elements
    [p, f, g, h, k, L] = elements.classical_to_equinoctial(a, e, incl, RAAN, omega, nu);
    x0_equinoctial = [p, f, g, h, k, L];
    [t_equinoctial_e, sol_equinoctial_e] = ode45(@propagate.equinoctial_pert, t, x0_equinoctial, opts_z, param);
    sol_equinoctial_e = sol_equinoctial_e';
    
    % Preallocate ECI position and velocity vectors
    r_eci_equinoctial_e = zeros(length(x0_equinoctial)/2, length(sol_equinoctial_e));
    v_eci_equinoctial_e = zeros(length(x0_equinoctial)/2, length(sol_equinoctial_e));
    
    % Convert equinoctial orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_equinoctial_e)
        p_sol = sol_equinoctial_e(1, i);
        f_sol = sol_equinoctial_e(2, i);
        g_sol = sol_equinoctial_e(3, i);
        h_sol = sol_equinoctial_e(4, i);
        k_sol = sol_equinoctial_e(5, i);
        L_sol = sol_equinoctial_e(6, i);
        [r_eci_equinoctial_e(:, i), v_eci_equinoctial_e(:, i)] = elements.equinoctial_to_cartesian(p_sol, f_sol, g_sol, h_sol, k_sol, L_sol, param.mu_earth);
    end
    sol_eci_equinoctial_e(:, :, z) = [r_eci_equinoctial_e; v_eci_equinoctial_e]; % Combine cartesian state vector in ECI
    
    % Simulate trajectory using milankovitch orbit elements
    [h_vec, e_vec, L] = elements.classical_to_milankovitch(a, e, incl, RAAN, omega, nu, param.mu_earth);
    x0_milankovitch_e = [h_vec; e_vec; L];
    [t_milankovitch_e, sol_milankovitch_e] = ode45(@propagate.milankovitch_pert, t, x0_milankovitch_e, opts_z, param);
    sol_milankovitch_e = sol_milankovitch_e';
    
    h_sol_milankovitch_e = sol_milankovitch_e(1:3, :); % Simulate milankovitch angular momentum vector (km^2/s)
    e_sol_milankovitch_e = sol_milankovitch_e(4:6, :); % Simulate milankovitch eccentricity vector
    
    milankovitch_constraint_e = zeros(1, length(sol_milankovitch_e), length(opts)); % Preallocate constraint variable
    
    % Calculate milankovitch constraint over time
    for i = 1:length(sol_milankovitch_e)
        milankovitch_constraint_e(1, i, z) = dot(h_sol_milankovitch_e(:, i), e_sol_milankovitch_e(:, i));
    end
    
    % Preallocate ECI position and velocity vectors
    r_eci_milankovitch_e = zeros(floor(length(x0_milankovitch_e)/2), length(sol_milankovitch_e));
    v_eci_milankovitch_e = zeros(floor(length(x0_milankovitch_e)/2), length(sol_milankovitch_e));
    
    % Convert equinoctial orbit elements from simulation to cartesian states in
    % ECI
    for i = 1:length(sol_milankovitch_e)
        h_sol_vec = sol_milankovitch_e(1:3, i);
        e_sol_vec = sol_milankovitch_e(4:6, i);
        L_sol = sol_milankovitch_e(7, i);
        [r_eci_milankovitch_e(:, i), v_eci_milankovitch_e(:, i)] = elements.milankovitch_to_cartesian(h_sol_vec, e_sol_vec, L_sol, param.mu_earth);
    end
    sol_eci_milankovitch_e(:, :, z) = [r_eci_milankovitch_e; v_eci_milankovitch_e]; % Combine cartesian state vector in ECI
end

% Plot orbit simulations on one plot

for z = 1:length(opts)
    figure(z + 3*length(opts) + 3)
    plot3(sol_cart_d(:, 1, z), sol_cart_d(:, 2, z), sol_cart_d(:, 3, z))
    hold on
    plot3(sol_eci_classical_e(1, :, z), sol_eci_classical_e(2, :, z), sol_eci_classical_e(3, :, z))
    plot3(sol_eci_equinoctial_e(1, :, z), sol_eci_equinoctial_e(2, :, z), sol_eci_equinoctial_e(3, :, z))
    plot3(sol_eci_milankovitch_e(1, :, z), sol_eci_milankovitch_e(2, :, z), sol_eci_milankovitch_e(3, :, z))
    title('Orbit Simulation Under J2, SRP, and Drag Perturbations - Sean Bohne', 'FontSize', 8)
    xlabel('x Position (km)')
    ylabel('y Position (km)')
    zlabel('z Position (km)')
    legend('Cartesian', 'Keplerian', 'Equinoctial', 'Milankovitch', 'location', 'best')
    axis equal
    grid on
end

% Plot milankovitch constraint over time
figure
for z = 1:length(opts)
    plot(t_milankovitch_e/3600, milankovitch_constraint_e(1, :, z))
    hold on
end
title('Milankovitch Constraint Over Time - Sean Bohne', 'FontSize', 8)
xlabel('Time (hrs)')
ylabel('Milankovitch Constraint')
grid on
legend('1e-12', '1e-9', '1e-6', '1e-3', 'location', 'southeast')