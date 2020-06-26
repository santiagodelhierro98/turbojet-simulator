function [Tt, Pt, ht, S_T, Eta_Overall, C_TS] = Turbojet_Sim_2(M0, T0, P0, Tt4, TauC, epsilon_i, epsilon_b, epsilon_n, eta_cp, eta_tp, fi, x, R, hf0)

% Data Arrays
Tt = [];
Pt = [];
ht = [];

Cp_air = cp_air(T0);
gamma = Cp_air/(Cp_air-R);

theta0 = (1+((gamma-1)/2)*M0^2);
Tt0 = T0*theta0;

phi_t0 = phi_air(Tt0);
phi_0 = phi_air(T0);
Pt0 = P0*exp(phi_t0 - phi_0);
ht0 = h_air(Tt0);

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
phi_t3 = phi_air(Tt3);
Pt3 = Pt2*exp((phi_t3 - phi_t2)*eta_cp);
ht3 = h_air(Tt3);

Tt = [Tt; Tt3];
Pt = [Pt; Pt3];
ht = [ht; ht3];

%% BURNER
Pt4 = Pt3*(1 - epsilon_b);
h4_air = h_air(Tt4);
ht4_fuel = h_fuel(Tt4);
phi_t4_air = phi_air(Tt4);
phi_t4_fuel = phi_fuel(Tt4);

h_fuel_T4 = h_fc(Tt4)*1000;
alpha = (h4_air - ht3)/h_fuel_T4;
alpha_ = alpha*(1 - x);

phi_t4 = (phi_t4_air + alpha*phi_t4_fuel)/(1 + alpha); 
ht4 = (h4_air + alpha*ht4_fuel)/(1 + alpha);

Tt = [Tt; Tt4];
Pt = [Pt; Pt4];
ht = [ht; ht4];

%% TURBINE
ht5 = ht4 - (ht3 - ht2)/((1 + alpha)*(1 - x));
ht5_mix = (ht5*(1 + alpha)*(1 - x) + ht3*x)/((1 + alpha)*(1 - x) + x);

function F = ht5_(x)
R = 287.15;
F(1) = R*(3.5*x - 1.4e-5*(x^2) + 7.467e-9*(x^3) + 3090/(exp(3090/x) - 1)) - ht5_mix;
end

x0 = Tt4;
Tt5 = fsolve(@ht5_, x0);

phi_t5 = (phi_air(Tt5) + alpha*phi_fuel(Tt5))/(1 + alpha);
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

pi_9 = Pt9/P0;
P9 = Pt9/pi_9;
phi_t9=phi_t5; % bc T9=T5
phi9i = phi_t9- exp(Pt9/P9);

function F = T9_(x)
F(1) = (3.5*log(x) - 2.8e-5*(x) + 1.12e-8*(x^2) + (3090/(x*(exp(3090/x) - 1))) - log((exp(3090/x) - 1)/exp(3090/x))) - phi_t9;
end
x0 = Tt9;
T9_is = fsolve(@T9_, x0);

cp_n = (cp_air(Tt9) + alpha*cp_fuel(Tt9))/(1+alpha);
V9 = fi*sqrt(2*cp_n*Tt9*((1-pi_9^((1-gamma)/gamma))));
v9 = V9/sqrt(2*Cp_air*T0);
V0 = M0*sqrt(gamma*R*T0);
v0 = V0/sqrt(2+Cp_air*T0);

S_T = (1 + alpha_)*V9 - V0;
Eta_Overall = S_T*V0/(alpha_*h_fuel_T4);
C_TS = (alpha_/S_T)*1e5; %% g/kN*s
Eta_Propulsive = 2*v0/(v9+v0);
Eta_Thermal = (V9^2 - V0^2)/(2*alpha_*hf0);
end