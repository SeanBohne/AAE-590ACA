function [r_eci, v_eci] = equinoctial_to_cartesian(p, f, g, h, k, L, mu)
% This function converts equinoctial orbit elements to cartesian states in
% ECI

% Calculate parameters for conversion
alpha_squared = h^2 - k^2;
s_squared = 1 + h^2 + k^2;
w = 1 + f*cos(L) + g*sin(L);
r = p/w;

r_eci = (r/s_squared)*[cos(L) + alpha_squared*cos(L) + 2*h*k*sin(L);
                       sin(L) - alpha_squared*sin(L) + 2*h*k*cos(L);
                       2*(h*sin(L) - k*cos(L))]; % ECI position vector (km)

v_eci = (1/s_squared*sqrt(mu/p))*[-(sin(L) + alpha_squared*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + alpha_squared*g);
                                  -(-cos(L) + alpha_squared*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + alpha_squared*f);
                                  2*(h*cos(L) + k*sin(L) + f*h + g*k)]; % ECI velocity vectory (km/s)
end