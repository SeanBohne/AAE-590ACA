function dxdt = classical(t, x, mu)
% This function propogates classical orbit elements using unperturbed 
% two-body motion

dxdt = zeros(6, 1);

n = sqrt(mu/x(1)^3); % Mean motion (rad/s)

dxdt(6) = n;
end