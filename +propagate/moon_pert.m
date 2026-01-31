function [a_moon, a_moon_eci] = moon_pert(t, r_vec, v_vec, GM_star, l_star, mu_grav)
% This function calculates the third body perturbation acceleration vector 
% in ECI and RTN frames

w_moon = sqrt(GM_star/l_star^3); % Moon angular velocity (rad/s)
r_moon = l_star; % Earth-Moon distance (km)

theta = w_moon*t; % Calculate angular position (rad)
r_moon_hat = [cos(theta);
              sin(theta);
              0]; % Moon position unit vector in ECI
r_vec_moon = r_moon*r_moon_hat; % Moon position vector in ECI (km)

r_rel_moon = r_vec - r_vec_moon; % Position vector between spacecraft and third body (km)
r_rel = norm(r_rel_moon); % Position magnitude between spacecraft and third body (km)

R = rotations.ECI_to_RTN(r_vec, v_vec); % Rotation matrix from ECI to RTN

a_moon_eci = -mu_grav*(r_rel_moon/r_rel^3 + r_vec_moon/r_moon^3); % Moon perturbation in ECI (km/s^2)
a_moon = R*a_moon_eci; % Moon perturbation in RTN (km/s^2)
end