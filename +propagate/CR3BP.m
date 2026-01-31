function dxdt = CR3BP(t, x, param)
% This function propogates cartesian states using unperturbed motion in
% CR3BP

% Store state to meaningful variables
x_pos = x(1);
y_pos = x(2);
z_pos = x(3);
vx = x(4);
vy = x(5);
vz = x(6);
mu_mass = param.mu_mass;

% Distance from spacecraft to each body
r13 = sqrt((x_pos + mu_mass)^2 + y_pos^2 + z_pos^2);
r23 = sqrt((x_pos - 1 + mu_mass)^2 + y_pos^2 + z_pos^2);

dxdt = zeros(6, 1);

dxdt(1) = vx;
dxdt(2) = vy;
dxdt(3) = vz;
dxdt(4) = 2*vy + x_pos - (1 - mu_mass)*(x_pos + mu_mass)/r13^3 - mu_mass*(x_pos - 1 + mu_mass)/r23^3;
dxdt(5) = -2*vx + y_pos - (1 - mu_mass)*y_pos/r13^3 - mu_mass*y_pos/r23^3;
dxdt(6) = -(1 - mu_mass)*z_pos/r13^3 - mu_mass*z_pos/r23^3;
end