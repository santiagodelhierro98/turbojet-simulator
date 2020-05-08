%% Turbojet Simulator
clear all;
close all;

%% Inputs
M0 = 2; % Mach 0
H = 11e3; % Altitude in meters
T0 = 216.5; % Temperature 0 in Kelvin
P0 = 22600; % P0 at 11km according to ISA
Tt4 = 1373; % Total Temperature 4 in Kelvin
TauC = 2; % Compressor T ratio Tt3/Tt2
epsilon_i = 0.075; epsilon_b = 0.06; epsilon_n = 0.02; % Pressure loss coefficients
eta_cp = 0.88; eta_tp = 0.93; % Politropic efficiencies
fi = 0.98; % Nozzle velocity coefficient
x = 0.07; % Turnine cooling bleed
R = 287.15; % Gas cte in J/kgK
hf0 = 4.3095e7; % Fuel lower heating value J/kg

%% Simulation

A = 3090/T0;
a=2.8e-5;
Cp_air = R*(3.5 - (2.8e-5)*T0 + (2.24e-8)*T0^2 + (A^2)*(exp(A)/(exp(A)-1)^2));
gamma = Cp_air/(Cp_air-R);

theta_0 = (1+((gamma-1)/2)*M0^2); % Tt0/T0
delta_0 = theta_0^(gamma/(gamma-1)); % Pt0/P0
h0 = R*(3.5*T0 - (T0^2)*1.4e-5 + (T0^3)*7.467e-9 + 3090/(A-1));
V0=M0*sqrt(gamma*R*T0);

% Intake: no work or heat input

theta_2 = theta_0; % Tt2/T0
delta_2 = (1-epsilon_i)*delta_0; % Pt2/P0
Pt2 = delta_0*P0;
h2 = h0;

% Compressor:

theta_3 = TauC*theta_0; % Tt3/T0
pi_c = TauC^(eta_cp*(gamma/(gamma-1))); % Pt3/Pt2
delta_3 = pi_c*delta_2; % Pt3/P0
Pt3=delta_3*P0;
Tt3 = theta_3*T0;
A = 3090/Tt3;
h3 = R*(3.5*Tt3 - (Tt3^2)*1.4e-5 + (Tt3^3)*7.467e-9 + 3090/(A-1));

% Burner:

pi_b = 1-epsilon_b; % Pt4/Pt3
theta_4 = Tt4/T0;
delta_4 = pi_b*delta_3; % Pt4/P0
Pt4 = delta_4*P0;
A = 3090/Tt4;
h4_air = R*(3.5*Tt4 - (Tt4^2)*1.4e-5 + (Tt4^3)*7.467e-9 + 3090/(A-1));
d_hfc = R*(-1607.2 + 4.47659*Tt4 + 4.00997e-3*(Tt4^2) - 6.12432e-7*(Tt4^3));
h4_fuel = hf0 - (d_hfc);

A = 3090/Tt4;
Cp4_air = R*(3.5 - (2.8e-5)*Tt4 + (2.24e-8)*Tt4^2 + (A^2)*(exp(A)/(exp(A)-1)^2));
Cp4_fuel = R*(4.47659 + 8.01994e-3*Tt4 - 1.873e-6*(Tt4^2));

alpha = (h4_air - h3)/h4_fuel;
%alpha = 0.5*((Cp4_fuel-1)+sqrt((1-Cp4_fuel)^2 + 4*Cp4_air*((Tt4 - Tt3)/h4_fuel)));
alpha_ = alpha*(1-x); %effective richness

h4 = (h4_air+alpha*h4_fuel)/(R+R*alpha);

% Bleeding

h5_bleed = h4 - (h3 - h2)/((1 - alpha)*(1 - x));

% Turbine:

h5 = (h5_bleed*(1 + alpha)*(1 - x) +h3*x)/((1 + alpha)*(1 - x) + x);
TauT = 1 - (theta_0/theta_4)*(TauC - 1);
theta_5 = TauT*theta_4; % Tt5/T0
pi_t = TauT^(gamma/(eta_tp*(gamma - 1)));
delat_5 = pi_t*delta_4; % Pt5/P0

% Nozzle:

theta_9 = theta_5; % Tt9/T0
delta_9 = (1 - epsilon_i)*(1 - epsilon_b)*(1 - epsilon_n)*theta_0^(gamma/(gamma - 1))*TauC^(gamma*eta_cp/(gamma - 1))*(1 - (theta_0/theta_4)*(TauC - 1))^(gamma/(eta_tp*(gamma - 1))); % Pt9/P0
epsilon_T = ((gamma - 1)/gamma)*(epsilon_i + epsilon_b + epsilon_n);
v9 = fi*sqrt(theta_4*(1-(1+epsilon_T)/((theta_0*TauC^eta_cp)*(1 - (theta_4/theta_0)*(TauC - 1))^((1 - eta_tp)/eta_tp))) - theta_0*(TauC - 1));
v0 = sqrt(theta_0 - 1);
V9 = sqrt(2*Cp_air*T0)*sqrt(theta_4 - (TauC*theta_0 - theta_0) - 1);

% Performances:

S_T = V9*(1+alpha_) - V0;
Eta_Overall = (S_T*V0)/(alpha_*h4_fuel);
C_TS = alpha_/S_T;

fprintf("Specific Thrust = " + S_T + "\n");
fprintf("Overall Efficiency = " + Eta_Overall + "\n");
fprintf("Specific Fuel Consumption = " + C_TS + "\n");