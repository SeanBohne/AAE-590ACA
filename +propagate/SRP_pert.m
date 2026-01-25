function [a_SPR, a_SRP_eci] = SRP_pert(t, r_vec, v_vec, A_m, G0, omega_sun, AU)
% This function calculates the canon-ball solar radiation pressure 
% perturbation acceleration vector in ECI and RTN frames

theta = omega_sun*t;
s_hat = [cos(theta);
           sin(theta);
           0]; % Sun unit vector

R = rotations.ECI_to_RTN(r_vec, v_vec); % Rotation matrix from ECI to RTN

a_SRP_eci = -A_m*G0/AU^2*s_hat; % SRP perturbation vector in ECI (km/s^2)
a_SPR = R*a_SRP_eci; % SRP perturbation in RTN (km/s^2)
end