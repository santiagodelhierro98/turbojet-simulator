function [Tt, Pt, ht, S_T, Eta_Overall, C_TS] = Turbojet_Sim_2(M0, T0, P0, Tt4, TauC, epsilon_i, epsilon_b, epsilon_n, eta_cp, eta_tp, fi, x, R, hf0)

% Data Arrays
Tt = [];
Pt = [];
ht = [];

A = 3090/T0;
Cp_air = R*(3.5 - (2.8e-5)*T0 + (2.24e-8)*T0^2 + (A^2)*(exp(A)/(exp(A)-1)^2));
gamma = Cp_air/(Cp_air-R);

theta0 = (1+((gamma-1)/2)*M0^2);
Tt0 = T0*theta0;
A = exp(3090/Tt0);
phi_0 = 3.5*log(Tt0) - 2.8e-5*Tt0 + 1.12e-8*(Tt0^2) + 3090/(Tt0*(A - 1)) - log((A - 1)/A);
Pr0 = exp(phi_0);
Pt0 = P0*theta0^(gamma/(gamma - 1));
ht0 = R*(3.5*T0 - (T0^2)*1.4e-5 + (T0^3)*7.467e-9 + 3090/(A-1));

Tt = [Tt; Tt0];
Pt = [Pt; Pt0];
ht = [ht; ht0];

% INTAKE
Tt2 = Tt0;
A = exp(3090/Tt2);
phi_2 = 3.5*log(Tt2) - 2.8e-5*Tt2 + 1.12e-8*(Tt2^2) + 3090/(Tt2*(A - 1)) - log((A - 1)/A);
Pr2 = Pr0;
Pt2 = Pt0*(1-epsilon_i);
ht2 = ht0;

Tt = [Tt; Tt2];
Pt = [Pt; Pt2];
ht = [ht; ht2];

% COMPRESSOR
Tt3 = TauC*Tt2;
A = exp(3090/Tt3);
phi_3 = 3.5*log(Tt3) - 2.8e-5*Tt3 + 1.12e-8*(Tt3^2) + 3090/(Tt3*(A - 1)) - log((A - 1)/A);
Pr3 = Pr2*exp(phi_3 - phi_2);
Pt3 = Pt2*(Pr3/Pr2)^eta_cp;
ht3 = R*(3.5*Tt3 - (Tt3^2)*1.4e-5 + (Tt3^3)*7.467e-9 + 3090/(A-1));

Tt = [Tt; Tt3];
Pt = [Pt; Pt3];
ht = [ht; ht3];

% BURNER
Pt4 = Pt3*(1 - epsilon_b);
A = exp(3090/Tt4);
h4_air = R*(3.5*Tt4 - (Tt4^2)*1.4e-5 + (Tt4^3)*7.467e-9 + 3090/(A-1));
h4_fuel = R*(-149.054 + 4.47659*Tt4 - 4.00997e-3*(Tt4^2) - 6.12432e-7*(Tt4^3));
phi_4_air = 3.5*log(Tt4) - 2.8e-5*Tt4 + 1.12e-8*(Tt4^2) + 3090/(Tt4*(A - 1)) - log((A - 1)/A);
phi_4_fuel = 4.47659*log(Tt4) + 8.01994e-3*Tt4 + 9.19648e-7*(Tt4^2);

h_fuel_T4 = hf0 - (-1607.2 + 4.47659*Tt4 + 4.00997e-3*(Tt4^2) - 6.12432e-7*(Tt4^3));
alpha = (h4_air - ht3)/h_fuel_T4;
alpha_ = alpha*(1 - x);

phi_4 = (phi_4_air + alpha*phi_4_fuel)/(1 + alpha); 
Pr4 = Pr3*exp(phi_4 - phi_3);
ht4 = (h4_air + alpha*h4_fuel)/(1 + alpha);

Tt = [Tt; Tt4];
Pt = [Pt; Pt4];
ht = [ht; ht4];

% TURBINE
ht5 = ht4 - (ht3 - ht2)/((1 + alpha)*(1 - x));
ht5_mix = (ht5*(1 + alpha)*(1 - x) + ht3*x)/((1 + alpha)*(1 - x) + x);

fun = R*(3.5*x - 1.4e-5*(x^2) + 7.467e-9*(x^3) + 3090/(exp(3090/x) - 1));
Tt5 = fsolve(fun, 983.28768, ht5);

end