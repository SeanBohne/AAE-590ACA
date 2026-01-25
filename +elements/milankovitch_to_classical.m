function [a, e, incl, RAAN, omega, nu, mu] = milankovitch_to_classical(h_vec, e_vec, L, mu)
% This function converts milankovitch orbit elements to classical orbit
% elements

h = norm(h_vec); % Specific angular momentum vector magnitude (km^2/s)
e = norm(e_vec); % Eccentricity

p = h^2/mu; % Semilatus rectum
a = p/(1 - e^2); % Semi-major axis (km)
incl = acos(h_vec(3)/h); % Inclination (rad)
RAAN = atan2(h_vec(1), -h_vec(2)); % Right ascension of ascending node (rad)
RAAN_omega = atan2(e_vec(2), e_vec(1)); % Longitude of periapsis (rad)
omega = RAAN_omega - RAAN; % Argument of periapsis (rad)
nu = L - RAAN_omega; % True anomaly (rad)
end