function dxdt = equinoctial(t, x, mu)
% This function propogates equinoctial orbit elements using unperturbed 
% two-body motion

dxdt = zeros(6, 1);
q = 1 + x(2)*cos(x(6)) + x(3)*sin(x(6)); % ODE parameter

dxdt(6) = sqrt(mu*x(1))*(q/x(1))^2;
end