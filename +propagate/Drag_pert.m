function [a_Drag, a_Drag_eci] = Drag_pert(r_vec, v_vec, rho_ref, dens_scale, h0, CD, A_m, w_earth)
% This function calculates the drag perturbation acceleration vector in ECI
% and RTN frames

% Calculate radial distance and relative velocity to atmosphere
r = norm(r_vec);
v_rel = v_vec - cross([0; 0; w_earth], r_vec);
v_mag = norm(v_rel);

R = rotations.ECI_to_RTN(r_vec, v_vec); % Rotation matrix from ECI to RTN

dens = rho_ref*exp(-dens_scale*(r - h0)); % Exponential density model (kg/km^3)
a_Drag_eci = -0.5*dens*CD*A_m*v_mag*v_rel; % Drag perturbation in ECI (km/s^2)
a_Drag = R*a_Drag_eci; % Drag perturbation in RTN (km/s^2)
end