function [r_CR3BP, v_CR3BP] = ECI_to_CR3BP(t, r_vec, v_vec, param)
% This function converts position and velocity vectors from dimensional 
% ECI to dimensional CR3BP

C = @(a) cos(a);
S = @(a) sin(a);

% Extract parameter values
l_star = param.l_star;
GM_star = param.GM_star;
t_star = sqrt(l_star^3/GM_star);
mu_mass = param.mu_mass;
w_moon = 1/t_star;
w_moon_vec = [0; 0; w_moon];

theta = w_moon*t; % Calculate angular position (rad)
R = [C(theta),  S(theta),  0;
     -S(theta), C(theta),  0;
     0,         0,         1]; % Rotation matrix

% Transform position and velocity vectors
r_CR3BP_rot = R*(r_vec);
v_CR3BP_rot = R*v_vec;

% Translate position and velocity vectors
r_CR3BP = r_CR3BP_rot + [-mu_mass*l_star; 0; 0]; 
v_CR3BP = v_CR3BP_rot - cross(w_moon_vec, r_CR3BP_rot);
end