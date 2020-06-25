function [cp] = cp_fuel(T)

cp = 287.15*(4.47659 + (8.01994e-3)*T - (1.8373e-6)*T*T);
end