function [a, e, incl, RAAN, omega, nu] = cartesian_to_classical(r_vec, v_vec, mu)
% This function converts cartesian states in ECI to classical orbit
% elements. Outputs angles in degrees.

r = norm(r_vec); % Distance magnitude (km)
r_hat = r_vec/r; % Distance unit vector
v_r = dot(v_vec, r_hat); % Radial velocity (km/s)

h_vec = cross(r_vec, v_vec); % Specific angular momentum vector (km^2/s)
h = norm(h_vec); % Specific angular momentum magnitude (km^2/s)
e_vec = cross(v_vec, h_vec)/mu - r_vec/r; % Eccentricity vector
e = norm(e_vec); % Eccentricity
p = h^2/mu; % Semilatus rectum

n_vec = cross([0, 0, 1], h_vec);
n = norm(n_vec);

a = p/(1 - e^2); % Semi-major axis (km)
incl = acos(h_vec(3) / h); % Inclination (rad)

if n_vec(2) >= 0
    RAAN = acos(n_vec(1)/n); % Right Ascension of ascending node (rad)
else
    RAAN = 360 - acos(n_vec(1)/n); % Right Ascension of ascending node (rad)
end

if e_vec(3) >= 0
    omega = acos(dot(e_vec, n_vec)/(e*n)); % Argument of periapsis (rad)
else
    omega = 360 - acos(dot(e_vec, n_vec)/(e*n)); % Argument of periapsis (rad)
end

if v_r >= 0
    nu = acos(dot(e_vec, r_vec)/(e*r)); % True anomaly (rad)
else
    nu = 360 - acos(dot(e_vec, r_vec)/(e*r)); % True anomaly (rad)
end
end