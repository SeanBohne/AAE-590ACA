function [a, e, incl, RAAN, omega, nu] = equinoctial_to_classical(p, f, g, h, k, L)
% This function converts equinoctial orbit elements to classical orbit
% elements

a = p/(1 - f^2 - g^2); % Semi-major axis (km)
e = sqrt(f^2 + g^2); % Eccentricity
incl = 2*atan(sqrt(h^2 + k^2)); % Inclination (rad)
RAAN = atan2(k, h); % Right ascension of ascending node (rad)
omega = atan2(g, f) - atan2(k, h); % Argument of periapsis (rad)
nu = L - (RAAN + omega); % True anomaly (rad)
end