function R = DCM_313(RAAN, incl, omega)
% This function calculates the rotation matrix based on a 3-1-3 rotation
% sequence. Takes degrees as input.

C = @(a) cos(a);
S = @(a) sin(a);

R = [C(RAAN) * C(omega) - S(RAAN) * S(omega) * C(incl), -C(RAAN) * S(omega) - S(RAAN) * C(omega) * C(incl), S(RAAN) * S(incl);
     S(RAAN) * C(omega) + C(RAAN) * S(omega) * C(incl), -S(RAAN) * S(omega) + C(RAAN) * C(omega) * C(incl), -C(RAAN) * S(incl);
     S(omega) * S(incl), C(omega) * S(incl), C(incl)]; % 3-1-3 rotation matrix
end