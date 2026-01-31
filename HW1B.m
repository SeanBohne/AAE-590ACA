clc; clear; close all;

%% Part A
% Initialize simulation parameters
param.l_star = 3.8475e5; % Earth-Moon distance (km)
param.GM_star = 4.035e5; % Earth-Moon barycenter GM (km^3/s^2)
param.mu_mass = 1.2151e-2; % Mass ratio
t_star = sqrt(param.l_star^3/param.GM_star); % Non-dimensionalizing time unit (s)
t_points = 15000; % Number of simulated time points

% Initial conditions
IC(1).x0_CR3BP = [1.2; 0; 0; 0; -1.06110124; 0];
IC(1).t_CR3BP = 6.20628;

IC(2).x0_CR3BP = [0.85; 0; 0.17546505; 0; 0.2628980369; 0];
IC(2).t_CR3BP = 2.5543991;

IC(3).x0_CR3BP = [0.05; -0.05; 0; 4; 2.6; 0];
IC(3).t_CR3BP = 15;

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); % ODE solver tolerance

% Preallocate state variables
x_CR3BP_dim = zeros(6, t_points, length(IC));
r_eci_a = zeros(3, t_points, length(IC));
v_eci_a = zeros(3, t_points, length(IC));

% Simulate CR3BP for each initial conditions and plot in dimensional system
for i = 1:length(IC)
    % ODE45 call for CR3BP simulation
    t_span_a = linspace(0, IC(i).t_CR3BP, t_points);
    [t_CR3BP, x_CR3BP] = ode45(@propagate.CR3BP, t_span_a, IC(i).x0_CR3BP, opts, param);
    x_CR3BP = x_CR3BP';
    
    % Dimensional quantities
    r_dim = param.l_star*x_CR3BP(1:3, :);
    v_dim = (param.l_star/t_star)*x_CR3BP(4:6, :);
    x_CR3BP_dim(:, :, i) = [r_dim; v_dim]; % Combine dimensional position and velocity (km, km/s)

    figure
    subplot(1, 2, 1)
    plot3(x_CR3BP_dim(1, :, i), x_CR3BP_dim(2, :, i), x_CR3BP_dim(3, :, i))
    hold on
    plot3(-param.l_star*param.mu_mass, 0, 0, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b') % Plot Earth location
    plot3(param.l_star*(1 - param.mu_mass), 0, 0, 'go', 'MarkerSize', 3, 'MarkerFaceColor', 'g') % Plot Moon location
    title('CR3BP Dimensional Trajectory in Synodic Frame', 'FontSize', 8)
    xlabel('x Position (km)')
    ylabel('y Position (km)')
    zlabel('z Position (km)')
    legend('Trajectory', 'Earth', 'Moon')
    grid on
    axis equal

    for j = 1:t_points
        t_eci_a = t_CR3BP(j);
        [r_eci_a(:, j, i), v_eci_a(:, j, i)] = rotations.CR3BP_to_ECI(t_eci_a, x_CR3BP(1:3, j), x_CR3BP(4:6, j), param);
    end

    subplot(1, 2, 2)
    plot3(r_eci_a(1, :, i), r_eci_a(2, :, i), r_eci_a(3, :, i))
    hold on
    plot3(0, 0, 0, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b')
    title('CR3BP Dimensional Trajectory in ECI Frame', 'FontSize', 8)
    xlabel('x Position (km)')
    ylabel('y Position (km)')
    zlabel('z Position (km)')
    % Plot Moon position at each time step
    r_moon_eci_a = zeros(3, t_points);
    for j = 1:t_points
        theta_moon = t_CR3BP(j);
        r_moon_eci_a(:, j) = param.l_star * [cos(theta_moon); sin(theta_moon); 0];
    end
    plot3(r_moon_eci_a(1, :), r_moon_eci_a(2, :), r_moon_eci_a(3, :), 'k')
    legend('Trajectory', 'Earth', 'Moon Path')
    axis equal
    sgtitle(['CR3BP Dimensional Trajectory in Synodic and ECI Frame - IC Set #', num2str(i)])
    grid on
end

%% Part B
param.use_3body = true; % Indicate that moon perturbations are used
param.use_J2 = false; % Indicate that J2 perturbations are not used
param.mu_grav = 4.9028e3; % Moon gravity parameter (km^3/s^2)
param.mu_grav_main = 3.986e5; % Earth gravity parameter (km^3/s^2)

% Preallocate state variables
x_3body_b = zeros(6, t_points, length(IC));
r_CR3BP_b = zeros(3, t_points, length(IC));
v_CR3BP_b = zeros(3, t_points, length(IC));
t_3body_b = zeros(t_points, length(IC));

% Simulate third body perturbations in ECI for each set of initial
% conditions
for i = 1:length(IC)
    % Convert CR3BP initial conditions to ECI
    [r_eci_b, v_eci_b] = rotations.CR3BP_to_ECI(0, IC(i).x0_CR3BP(1:3), IC(i).x0_CR3BP(4:6), param);
    x0_eci_b = [r_eci_b; v_eci_b];
    t_eci_b = t_star*IC(i).t_CR3BP;
    tspan_b = linspace(0, t_eci_b, t_points);

    % ODE45 call for third-body simulation
    [t_3body_intermediate, x_3body_intermediate] = ode45(@propagate.cartesian_pert, tspan_b, x0_eci_b, opts, param);
    x_3body_intermediate = x_3body_intermediate';
    x_3body_b(:, :, i) = x_3body_intermediate;
    t_3body_b(:, i) = t_3body_intermediate;

    figure
    subplot(1, 2, 2)
    plot3(x_3body_b(1, :, i), x_3body_b(2, :, i), x_3body_b(3, :, i))
    hold on
    plot3(0, 0, 0, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b')
    title('Trajectory With Third-Body Perturbations in ECI Frame', 'FontSize', 8)
    xlabel('x Position (km)')
    ylabel('y Position (km)')
    zlabel('z Position (km)')
    % Plot Moon position at each time step
    r_moon_eci_b = zeros(3, t_points);
    for j = 1:t_points
        theta_moon = t_3body_b(j, i)/t_star;
        r_moon_eci_b(:, j) = param.l_star * [cos(theta_moon); sin(theta_moon); 0];
    end
    plot3(r_moon_eci_b(1, :), r_moon_eci_b(2, :), r_moon_eci_b(3, :), 'k')
    legend('Trajectory', 'Earth', 'Moon Path')
    axis equal
    grid on

    for j = 1:t_points
        t_eci_b = t_3body_b(j, i);
        [r_CR3BP_b(:, j, i), v_CR3BP_b(:, j, i)] = rotations.ECI_to_CR3BP(t_eci_b, x_3body_b(1:3, j, i), x_3body_b(4:6, j, i), param);
    end

    subplot(1, 2, 1)
    plot3(r_CR3BP_b(1, :, i), r_CR3BP_b(2, :, i), r_CR3BP_b(3, :, i))
    hold on
    plot3(-param.l_star*param.mu_mass, 0, 0, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b') % Plot Earth location
    plot3(param.l_star*(1 - param.mu_mass), 0, 0, 'go', 'MarkerSize', 3, 'MarkerFaceColor', 'g') % Plot Moon location
    title('Trajectory With Third-Body Perturbations in Synodic Frame', 'FontSize', 8)
    xlabel('x Position (km)')
    ylabel('y Position (km)')
    zlabel('z Position (km)')
    legend('Trajectory', 'Earth', 'Moon')
    axis equal
    sgtitle(['Third-Body Perturbation Trajectory in Synodic and ECI Frame - IC Set #', num2str(i)])
    grid on
end

%% Part C
% Calculate differences in CR3BP and third body perturbation simulations
% and plot
for i = 1:length(IC)
    r_rot_diff_c = squeeze(x_CR3BP_dim(1:3, :, i) - r_CR3BP_b(1:3, :, i));
    r_eci_diff_c = squeeze(r_eci_a(1:3, :, i) - x_3body_b(1:3, :, i));

    v_rot_diff_c = squeeze(x_CR3BP_dim(4:6, :, i) - v_CR3BP_b(1:3, :, i));
    v_eci_diff_c = squeeze(v_eci_a(1:3, :, i) - x_3body_b(4:6, :, i));
    
    figure
    subplot(1, 3, 1)
    plot(t_3body_b(:, i)/(3600*24), r_rot_diff_c(1, :))
    title('x Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('x Position Difference (km)')
    grid on
    subplot(1, 3, 2)
    plot(t_3body_b(:, i)/(3600*24), r_rot_diff_c(2, :))
    title('y Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('y Position Difference (km)')
    grid on
    subplot(1, 3, 3)
    plot(t_3body_b(:, i)/(3600*24), r_rot_diff_c(3, :))
    title('z Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('z Position Difference (km)')
    sgtitle(['CR3BP and Third-Body Perturbation Differences in Synodic Frame - IC Set #', num2str(i)])
    grid on

    figure
    subplot(1, 3, 1)
    plot(t_3body_b(:, i)/(3600*24), v_rot_diff_c(1, :))
    title('x Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('x Velocity Difference (km/s)')
    grid on
    subplot(1, 3, 2)
    plot(t_3body_b(:, i)/(3600*24), v_rot_diff_c(2, :))
    title('y Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('y Velocity Difference (km/s)')
    grid on
    subplot(1, 3, 3)
    plot(t_3body_b(:, i)/(3600*24), v_rot_diff_c(3, :))
    title('z Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('z Velocity Difference (km/s)')
    sgtitle(['CR3BP and Third-Body Perturbation Differences in Synodic Frame - IC Set #', num2str(i)])
    grid on
    
    figure
    subplot(1, 3, 1)
    plot(t_3body_b(:, i)/(3600*24), r_eci_diff_c(1, :))
    title('x Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('x Position Difference (km)')
    grid on
    subplot(1, 3, 2)
    plot(t_3body_b(:, i)/(3600*24), r_eci_diff_c(2, :))
    title('y Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('y Position Difference (km)')
    grid on
    subplot(1, 3, 3)
    plot(t_3body_b(:, i)/(3600*24), r_eci_diff_c(3, :))
    title('z Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('z Position Difference (km)')
    sgtitle(['CR3BP and Third-Body Perturbation Differences in ECI Frame - IC Set #', num2str(i)])
    grid on

    figure
    subplot(1, 3, 1)
    plot(t_3body_b(:, i)/(3600*24), v_eci_diff_c(1, :))
    title('x Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('x Velocity Difference (km/s)')
    grid on
    subplot(1, 3, 2)
    plot(t_3body_b(:, i)/(3600*24), v_eci_diff_c(2, :))
    title('y Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('y Velocity Difference (km/s)')
    grid on
    subplot(1, 3, 3)
    plot(t_3body_b(:, i)/(3600*24), v_eci_diff_c(3, :))
    title('z Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('z Velocity Difference (km/s)')
    sgtitle(['CR3BP and Third-Body Perturbation Differences in ECI Frame - IC Set #', num2str(i)])
    grid on
end

%% Part D
param.use_3body = true; % Indicate that moon perturbations are used
param.use_J2 = true; % Indicate that J2 perturbations are used
param.r_earth = 6378.1; % Earth radius (km)
param.J2 = 1.0826e-3; % J2 perturbation parameter

% Preallocate state variables
x_3body_d = zeros(6, t_points, length(IC));
r_CR3BP_d = zeros(3, t_points, length(IC));
v_CR3BP_d = zeros(3, t_points, length(IC));
t_3body_d = zeros(t_points, length(IC));

% Simulate third body and J2 perturbations in ECI for each set of initial
% conditions
for i = 1:length(IC)
    % Convert CR3BP initial conditions to ECI
    [r_eci_d, v_eci_d] = rotations.CR3BP_to_ECI(0, IC(i).x0_CR3BP(1:3), IC(i).x0_CR3BP(4:6), param);
    x0_eci_d = [r_eci_d; v_eci_d];
    t_eci_d = t_star*IC(i).t_CR3BP;
    tspan_d = linspace(0, t_eci_d, t_points);

    % ODE45 call for third-body simulation
    [t_3body_intermediate, x_3body_intermediate] = ode45(@propagate.cartesian_pert, tspan_d, x0_eci_d, opts, param);
    x_3body_intermediate = x_3body_intermediate';
    x_3body_d(:, :, i) = x_3body_intermediate;
    t_3body_d(:, i) = t_3body_intermediate;

    figure
    subplot(1, 2, 2)
    plot3(x_3body_d(1, :, i), x_3body_d(2, :, i), x_3body_d(3, :, i))
    hold on
    plot3(0, 0, 0, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b')
    title('Trajectory With Third-Body and J2 Perturbations in ECI Frame', 'FontSize', 8)
    xlabel('x Position (km)')
    ylabel('y Position (km)')
    zlabel('z Position (km)')
    % Plot Moon position at each time step
    r_moon_eci_d = zeros(3, t_points);
    for j = 1:t_points
        theta_moon = t_3body_d(j, i)/t_star;
        r_moon_eci_d(:, j) = param.l_star * [cos(theta_moon); sin(theta_moon); 0];
    end
    plot3(r_moon_eci_d(1, :), r_moon_eci_d(2, :), r_moon_eci_d(3, :), 'k')
    legend('Trajectory', 'Earth', 'Moon Path')
    axis equal
    grid on

    for j = 1:t_points
        t_eci_d = t_3body_d(j, i);
        [r_CR3BP_d(:, j, i), v_CR3BP_d(:, j, i)] = rotations.ECI_to_CR3BP(t_eci_d, x_3body_d(1:3, j, i), x_3body_d(4:6, j, i), param);
    end

    subplot(1, 2, 1)
    plot3(r_CR3BP_d(1, :, i), r_CR3BP_d(2, :, i), r_CR3BP_d(3, :, i))
    hold on
    plot3(-param.l_star*param.mu_mass, 0, 0, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b') % Plot Earth location
    plot3(param.l_star*(1 - param.mu_mass), 0, 0, 'go', 'MarkerSize', 3, 'MarkerFaceColor', 'g') % Plot Moon location
    title('Trajectory With Third-Body and J2 Perturbations in Synodic Frame', 'FontSize', 8)
    xlabel('x Position (km)')
    ylabel('y Position (km)')
    zlabel('z Position (km)')
    legend('Trajectory', 'Earth', 'Moon')
    axis equal
    sgtitle(['Third-Body and J2 Perturbation Trajectory in Synodic and ECI Frame - IC Set #', num2str(i)])
    grid on
end

% Calculate differences in third-body perturbation simulations with and 
% without J2 perturbations and plot
for i = 1:length(IC)
    r_rot_diff_d = squeeze(r_CR3BP_b(1:3, :, i) - r_CR3BP_d(1:3, :, i));
    r_eci_diff_d = squeeze(x_3body_b(1:3, :, i) - x_3body_d(1:3, :, i));

    v_rot_diff_d = squeeze(v_CR3BP_b(1:3, :, i) - v_CR3BP_d(1:3, :, i));
    v_eci_diff_d = squeeze(x_3body_b(4:6, :, i) - x_3body_d(4:6, :, i));
    
    figure
    subplot(1, 3, 1)
    plot(t_3body_d(:, i)/(3600*24), r_rot_diff_d(1, :))
    title('x Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('x Position Difference (km)')
    grid on
    subplot(1, 3, 2)
    plot(t_3body_d(:, i)/(3600*24), r_rot_diff_d(2, :))
    title('y Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('y Position Difference (km)')
    grid on
    subplot(1, 3, 3)
    plot(t_3body_d(:, i)/(3600*24), r_rot_diff_d(3, :))
    title('z Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('z Position Difference (km)')
    sgtitle(['Third-Body and J2 Perturbation Differences in Synodic Frame - IC Set #', num2str(i)])
    grid on

    figure
    subplot(1, 3, 1)
    plot(t_3body_d(:, i)/(3600*24), v_rot_diff_d(1, :))
    title('x Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('x Velocity Difference (km/s)')
    grid on
    subplot(1, 3, 2)
    plot(t_3body_d(:, i)/(3600*24), v_rot_diff_d(2, :))
    title('y Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('y Velocity Difference (km/s)')
    grid on
    subplot(1, 3, 3)
    plot(t_3body_d(:, i)/(3600*24), v_rot_diff_d(3, :))
    title('z Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('z Velocity Difference (km/s)')
    sgtitle(['Third-Body and J2 Perturbation Differences in Synodic Frame - IC Set #', num2str(i)])
    grid on
    
    figure
    subplot(1, 3, 1)
    plot(t_3body_d(:, i)/(3600*24), r_eci_diff_d(1, :))
    title('x Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('x Position Difference (km)')
    grid on
    subplot(1, 3, 2)
    plot(t_3body_d(:, i)/(3600*24), r_eci_diff_d(2, :))
    title('y Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('y Position Difference (km)')
    grid on
    subplot(1, 3, 3)
    plot(t_3body_d(:, i)/(3600*24), r_eci_diff_d(3, :))
    title('z Position Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('z Position Difference (km)')
    sgtitle(['Third-Body and J2 Perturbation Differences in ECI Frame - IC Set #', num2str(i)])
    grid on

    figure
    subplot(1, 3, 1)
    plot(t_3body_d(:, i)/(3600*24), v_eci_diff_d(1, :))
    title('x Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('x Velocity Difference (km/s)')
    grid on
    subplot(1, 3, 2)
    plot(t_3body_d(:, i)/(3600*24), v_eci_diff_d(2, :))
    title('y Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('y Velocity Difference (km/s)')
    grid on
    subplot(1, 3, 3)
    plot(t_3body_d(:, i)/(3600*24), v_eci_diff_d(3, :))
    title('z Velocity Difference - Sean Bohne', 'FontSize', 8)
    xlabel('Time (days)')
    ylabel('z Velocity Difference (km/s)')
    sgtitle(['Third-Body and J2 Perturbation Differences in ECI Frame - IC Set #', num2str(i)])
    grid on
end