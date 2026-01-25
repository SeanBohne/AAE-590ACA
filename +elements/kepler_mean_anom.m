function M = kepler_mean_anom(nu, e)
% This function solves Kepler's equation to find mean anomaly given true
% anomaly

E = 2*atan2(sqrt(1 - e)*sin(nu/2), sqrt(1 + e)*cos(nu/2)); % Eccentric anomaly (rad)
M = E - e*sin(E); % Mean anomaly (rad)
end