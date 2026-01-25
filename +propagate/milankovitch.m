function dxdt = milankovitch(t, x, mu)
% This function propogates milankovitch orbit elements using unperturbed 
% two-body motion

dxdt = zeros(7, 1);

% Store states in meaningful variables
h_vec = x(1:3);
e_vec = x(4:6);
L = x(7);

h = norm(h_vec); % Angular momentum magnitude (km^2/s)

% Find radial distance (km)
[r_vec, ~] = elements.milankovitch_to_cartesian(h_vec, e_vec, L, mu);
r = norm(r_vec);

dxdt(7) = h/r^2;
end