function [r_eci, v_eci] = classical_to_cartesian(a, e, incl, RAAN, omega, nu, mu)
% This function converts classical orbit elements to cartesian states in
% ECI

r = a*(1 - e^2)/(1 + e*cos(nu)); % Radial position center of Earth (km)
r_pf = r* [cos(nu);
           sin(nu);
           0]; % Perifocal position (km)

h = sqrt(mu*a*(1 - e^2)); % Specific angular momentum (km^2/s)

v_pf = (mu/h)*[-sin(nu);
               e + cos(nu);
               0]; % Perifocal velocity (km/s)
R = rotations.DCM_313(RAAN, incl, omega); % Calculate rotation matrix to ECI

r_eci = R*r_pf; % Convert to ECI
v_eci = R*v_pf; % Convert to ECI
end