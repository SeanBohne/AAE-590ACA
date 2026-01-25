function [h_vec, e_vec, L, mu] = cartesian_to_milankovitch(r_vec, v_vec, mu)
% This function converts cartesian states in ECI to milankovitch orbit 
% elements

% Convert to classical orbit elements then to milankovitch orbit elements
[a, e, incl, RAAN, omega, nu] = cartesian_to_classical(r_vec, v_vec, mu);
[h_vec, e_vec, L] = classical_to_milankovitch(a, e, incl, RAAN, omega, nu, mu);
end