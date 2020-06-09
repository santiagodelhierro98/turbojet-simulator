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
phi_t0 = 3.5*log(Tt0) - 2.8e-5*Tt0 + 1.12e-8*(Tt0^2) + 3090/(Tt0*(A - 1)) - log((A - 1)/A);
A = exp(3090/T0);
phi_0 = 3.5*log(T0) - 2.8e-5*T0 + 1.12e-8*(T0^2) + 3090/(T0*(A - 1)) - log((A - 1)/A);
Pt0 = P0*exp(phi_t0 - phi_0);
ht0 = R*(3.5*Tt0 - (Tt0^2)*1.4e-5 + (Tt0^3)*7.467e-9 + 3090/(A-1));

Tt = [Tt; Tt0];
Pt = [Pt; Pt0];
ht = [ht; ht0];

%% INTAKE
Tt2 = Tt0;
phi_t2 = phi_t0;
Pt2 = Pt0*(1-epsilon_i);
ht2 = ht0;

Tt = [Tt; Tt2];
Pt = [Pt; Pt2];
ht = [ht; ht2];

%% COMPRESSOR
Tt3 = TauC*Tt2;
A = exp(3090/Tt3);
phi_t3 = 3.5*log(Tt3) - 2.8e-5*Tt3 + 1.12e-8*(Tt3^2) + 3090/(Tt3*(A - 1)) - log((A - 1)/A);
Pt3 = Pt2*exp((phi_t3 - phi_t2)*eta_cp);
ht3 = R*(3.5*Tt3 - (Tt3^2)*1.4e-5 + (Tt3^3)*7.467e-9 + 3090/(A - 1));

Tt = [Tt; Tt3];
Pt = [Pt; Pt3];
ht = [ht; ht3];

%% BURNER
Pt4 = Pt3*(1 - epsilon_b);
A = exp(3090/Tt4);
h4_air = R*(3.5*Tt4 - (Tt4^2)*1.4e-5 + (Tt4^3)*7.467e-9 + 3090/(A-1));
h4_fuel = R*(-149.054 + 4.47659*Tt4 + 4.00997e-3*(Tt4^2) - 6.12432e-7*(Tt4^3));
phi_t4_air = 3.5*log(Tt4) - 2.8e-5*Tt4 + 1.12e-8*(Tt4^2) + 3090/(Tt4*(A - 1)) - log((A - 1)/A);
phi_t4_fuel = 4.47659*log(Tt4) + 8.01994e-3*Tt4 + 9.19648e-7*(Tt4^2);

h_fuel_T4 = hf0 - (-1607.2 + 4.47659*Tt4 + 4.00997e-3*(Tt4^2) - 6.12432e-7*(Tt4^3));
alpha = (h4_air - ht3)/h_fuel_T4;
alpha_ = alpha*(1 - x);

phi_t4 = (phi_t4_air + alpha*phi_t4_fuel)/(1 + alpha); 
ht4 = (h4_air + alpha*h4_fuel)/(1 + alpha);

Tt = [Tt; Tt4];
Pt = [Pt; Pt4];
ht = [ht; ht4];

%% TURBINE
ht5 = ht4 - (ht3 - ht2)/((1 + alpha)*(1 - x));
ht5_mix = (ht5*(1 + alpha)*(1 - x) + ht3*x)/((1 + alpha)*(1 - x) + x);

x0 = Tt4;
Tt5 = fsolve(@ht5, x0);

A = exp(3090/Tt5);
phi_t5 = 3.5*log(Tt5) - 2.8e-5*Tt5 + 1.12e-8*(Tt5^2) + 3090/(Tt5*(A - 1)) - log((A - 1)/A);
Pt5 = Pt4*exp((phi_t5 - phi_t4)*(1/eta_tp));

Tt = [Tt; Tt5];
Pt = [Pt; Pt5];
ht = [ht; ht5_mix];


%% NOZZLE
Pt9 = Pt5*(1 - epsilon_n);
Tt9 = Tt5;
ht9 = ht5_mix;

Tt = [Tt; Tt9];
Pt = [Pt; Pt9];
ht = [ht; ht9];

phi_t9 = phi_t5;
pi_9 = Pt9/P0;
phi_9 = phi_t9 - pi_9;

x0 = Tt9;
T9_is = fsolve(@T9, x0);
A = exp(3090/T9_is);
h9_is = R*(3.5*T9_is - (T9_is^2)*1.4e-5 + (T9_is^3)*7.467e-9 + 3090/(A-1));

V9_is = sqrt(2*(ht9 - h9_is));
V9 = (V9_is*fi);
v9 = V9/sqrt(2*Cp_air*T0);
V0 = M0*sqrt(gamma*R*T0);
v0 = V0/sqrt(2+Cp_air*T0);

S_T = (1 + alpha_)*V9 - V0;
Eta_Overall = S_T*V0/(alpha_*h_fuel_T4);
C_TS = (alpha_/S_T)*1e5; %% g/kN*s
Eta_Propulsive = 2*v0/(v9+v0);
Eta_Thermal = (V9^2 - V0^2)/(2*alpha_*hf0);
end