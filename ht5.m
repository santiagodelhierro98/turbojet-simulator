function F = ht5(x)

R = 287.15;
ht5 = 8.290350423510503e+05;
F(1) = R*(3.5*x - 1.4e-5*(x^2) + 7.467e-9*(x^3) + 3090/(exp(3090/x) - 1)) - ht5;

end