function [phi] = phi_fuel(T)

A = exp(3090/T);
phi = 4.47659*log(T) + 8.01994e-3*T + 9.19648e-7*(T^2);

end