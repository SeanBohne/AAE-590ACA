function R = ECI_to_RTN(r_vec, v_vec)
% This function calculates the rotation matrix to convert states in ECI to
% RTN. Takes degrees as input.

% Unit vectors in RTN frame
r_hat = r_vec/norm(r_vec);
n_hat = cross(r_vec, v_vec)/norm(cross(r_vec, v_vec));
t_hat = cross(n_hat, r_hat);

R = [r_hat';
     t_hat';
     n_hat'];
end