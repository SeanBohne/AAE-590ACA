function [h_vec, e_vec, L] = classical_to_milankovitch(a, e, incl, RAAN, omega, nu, mu)
% This function converts classical orbit elements to milankovitch orbit
% elements

h = sqrt(mu*a*(1 - e^2)); % Specific angular momentum (km^2/s)
h_hat = [sin(RAAN)*sin(incl);
         -cos(RAAN)*sin(incl);
         cos(incl)]; % Specific angular momentum unit vector (km^2/s)

e_hat = [cos(omega)*cos(RAAN) - cos(incl)*sin(omega)*sin(RAAN);
         cos(omega)*sin(RAAN) + cos(incl)*sin(omega)*cos(RAAN);
         sin(omega)*sin(incl)]; % Eccentricity unit vector

h_vec = h*h_hat; % Specific angular momentum vector (km^2/s)
e_vec = e*e_hat; % Eccentricity vector
RAAN_omega = atan2(e_vec(2), e_vec(1)); % Longitude of periapsis (rad)
L = RAAN_omega + nu; % or project r_vec in orbital plane
end