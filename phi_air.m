function [phi] = phi_air(T)

% A = exp(3090/T);
% 
% a1 = 3.5*log(T);
% b1 = (2.8e-5)*T;
% c1 = (1.12e-8)*(T^2);
% d1 = 3090/(T*(A - 1));
% f1 = log((A - 1)/A);
% 
% phi = a1 - b1 + c1 + d1 - f1;

a2 = 3.5*log(T);
b2 = (2.8e-5)*T;
c2 = (1.12e-8)*(T^2);
d2 = 3090/(T*(exp(3090/T)-1));
f2 = log((exp(3090/T)-1)/(exp(3090/T)));

phi = a2 - b2 + c2 + d2 - f2;

end