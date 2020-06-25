function [h] = h_air(T)

A = exp(3090/T);
h = 287.15*(3.5*T - (T^2)*1.4e-5 + (T^3)*7.467e-9 + 3090/(A-1));

end