function [p, f, g, h, k, L] = cartesian_to_equinoctial(r_vec, v_vec, mu)
% This function converts cartesian states in ECI to equinoctial orbit 
% elements

% Convert to classial elements then to equinoctial
[a, e, incl, RAAN, omega, nu] = cartesian_to_classical(r_vec, v_vec, mu);
[p, f, g, h, k, L] = classical_to_equinoctial(a, e, incl, RAAN, omega, nu, mu);
end