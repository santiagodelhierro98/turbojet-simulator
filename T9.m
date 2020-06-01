function F = T9(x)

phi_9 = 17.172624259376637;
F(1) = (3.7*log(x) - 2.8e-5*x + 1.12e-8*(x^2) + 3090/(x*(exp(3090/x) - 1)) - log((exp(3090/x) - 1)/exp(3090/x))) - phi_9;

end