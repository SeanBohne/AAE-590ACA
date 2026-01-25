function [r_vec, v_vec] = milankovitch_to_cartesian(h_vec, e_vec, L, mu)
% This function converts milankovitch orbit elements to cartesian states in
% ECI

% Calculate angular momentum and eccentricity magnitudes
h = norm(h_vec);
e = norm(e_vec);

% calculate orthonormal basis vectors
h_hat = h_vec/h;
e_hat = e_vec/e;
e_perp_hat = cross(h_hat, e_hat);
e_perp_hat = e_perp_hat/norm(e_perp_hat); % Ensure vector has norm 1

nu = mod(L - atan2(e_vec(2), e_vec(1)), 360); % True anomaly (rad)
p = h^2/mu; % Semilatus rectum
r = p/(1 + e*cos(nu)); % Radial position (km)

r_vec = r*(cos(nu)*e_hat + sin(nu)*e_perp_hat); % ECI position vector (km)
v_vec = sqrt(mu/p)*(-sin(nu)*e_hat + (e + cos(nu))*e_perp_hat); % ECI velocity vector (km/s)
end