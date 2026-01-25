function dxdt = cartesian(t, x, mu)
% This function propogates cartesian states using unperturbed two-body
% motion

r = norm(x(1:3)); % Distance between orbited body and satellite (km)

dxdt(1:3) = x(4:6);
dxdt(4:6) = -mu*x(1:3)/(r^3);
dxdt = dxdt';
end