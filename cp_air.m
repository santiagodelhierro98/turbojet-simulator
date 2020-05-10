function [cp] = cp_air(T)

r = 287.15;

cp_r = 3.5 - (2.8e-5)*T + (2.24e-8)*(T^2) + ((3090/T)^2)*(exp(3090/T)/(exp(3090/T)-1)^2);

cp = cp_r * r;

end
