function [r_eci, v_eci] = CR3BP_to_ECI(t, r_vec, v_vec, param)
% This function converts position and velocity vectors from dimensionless 
% CR3BP to dimensional ECI

C = @(a) cos(a);
S = @(a) sin(a);

% Extract parameter values
l_star = param.l_star;
GM_star = param.GM_star;
t_star = sqrt(l_star^3/GM_star);
mu_mass = param.mu_mass;
w_moon = 1;
w_moon_vec = [0; 0; w_moon];

theta = w_moon*t; % Calculate angular position (rad)
R = [C(theta), -S(theta), 0;
     S(theta), C(theta),  0;
     0,        0,         1]; % Rotation matrix

% Transform position and velocity vectors
r_rot = r_vec + [mu_mass; 0; 0];
r_eci_nondim = R*(r_rot);
v_eci_nondim = R*(v_vec + cross(w_moon_vec, r_rot));

% Dimensionalize position and velocity vectors
r_eci = r_eci_nondim*l_star; 
v_eci = v_eci_nondim*(l_star/t_star);
end