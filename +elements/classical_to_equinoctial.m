function [p, f, g, h, k, L] = classical_to_equinoctial(a, e, incl, RAAN, omega, nu)
% This function converts classical orbit elements to equinoctial orbit
% elements

p = a*(1 - e^2); % Semilatus rectum
f = e*cos(omega + RAAN);
g = e*sin(omega + RAAN);
h = tan(incl/2)*cos(RAAN);
k = tan(incl/2)*sin(RAAN);
L = RAAN + omega + nu; % True longtidue (rad)
end