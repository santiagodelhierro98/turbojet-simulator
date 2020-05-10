function [hf_Tt4] = h_fuel(T)

hf0 = 4.3095e7;
r = 287.15;
Delta_hf_entre_r = -1607.2 + 4.47659*T + 4.00997e-3*(T^2) - 6.12432e-7*(T)^3;
Delta_hf = Delta_hf_entre_r*r;
hf_Tt4 = hf0 - Delta_hf;

hf_Tt4 = hf_Tt4/1000;

end