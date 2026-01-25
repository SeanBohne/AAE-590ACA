function [a_J2, a_J2_eci] = J2_pert(r_vec, v_vec, mu, r0, J2)
% This function calculates the J2 perturbation acceleration vector in ECI 
% and RTN frames

r = norm(r_vec); % Distance between orbited body and satellite (km)
R = rotations.ECI_to_RTN(r_vec, v_vec); % Rotation matrix from ECI to RTN

% Precompute common terms
r2 = r^2;
z2 = r_vec(3)^2;
z2_r2 = z2/r2;
factor = -(3*mu*J2*r0^2)/(2*r^5);

a_J2_eci = factor*[(1 - 5*z2_r2)*r_vec(1);
               (1 - 5*z2_r2)*r_vec(2);
               (3 - 5*z2_r2)*r_vec(3)]; % J2 perturbation in ECI

a_J2 = R*a_J2_eci; % J2 perturbation in RTN
end