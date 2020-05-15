%% Turbojet Function

function [Tt, Pt, Ht, S_T, Eta_Overall, C_TS] = Turbojet_Sim(M0, T0, P0, Tt4, TauC, epsilon_i, epsilon_b, epsilon_n, eta_cp, eta_tp, fi, x, R, hf0)

% Data Arrays
Tt = [];
Pt = [];
Ht = [];

A = 3090/T0;
Cp_air = R*(3.5 - (2.8e-5)*T0 + (2.24e-8)*T0^2 + (A^2)*(exp(A)/(exp(A)-1)^2));
gamma = Cp_air/(Cp_air-R);

theta_0 = (1+((gamma-1)/2)*M0^2); % Tt0/T0
delta_0 = theta_0^(gamma/(gamma-1)); % Pt0/P0
h0 = R*(3.5*T0 - (T0^2)*1.4e-5 + (T0^3)*7.467e-9 + 3090/(A-1));
V0=M0*sqrt(gamma*R*T0);

Tt = [Tt; theta_0];
Pt = [Pt; delta_0];
Ht = [Ht; h0/h0];
% Intake: no work or heat input

theta_2 = theta_0; % Tt2/T0
delta_2 = (1-epsilon_i)*delta_0; % Pt2/P0
Pt2 = delta_2*P0;
h2 = h0;

Tt = [Tt; theta_2];
Pt = [Pt; delta_2];
Ht = [Ht; h2/h0];
% Compressor:

theta_3 = TauC*theta_0; % Tt3/T0
pi_c = TauC^(eta_cp*(gamma/(gamma-1))); % Pt3/Pt2
delta_3 = pi_c*delta_2; % Pt3/P0
Pt3 = delta_3*P0;
Tt3 = theta_3*T0;
B = 3090/Tt3;
h3 = R*(3.5*Tt3 - (Tt3^2)*1.4e-5 + (Tt3^3)*7.467e-9 + 3090/(B-1));

Tt = [Tt; theta_3];
Pt = [Pt; delta_3];
Ht = [Ht; h3/h0];
% Burner:

theta_4 = Tt4/T0;
pi_b = 1-epsilon_b; % Pt4/Pt3
delta_4 = pi_b*delta_3; % Pt4/P0
Pt4 = delta_4*P0;

C = 3090/Tt4;
h4_air = R*(3.5*Tt4 - (Tt4^2)*1.4e-5 + (Tt4^3)*7.467e-9 + 3090/(C-1));
d_hfc = R*(-1607.2 + 4.47659*Tt4 + 4.00997e-3*(Tt4^2) - 6.12432e-7*(Tt4^3));
h4_fuel = hf0 - (d_hfc);

Cp4_air = R*(3.5 - (2.8e-5)*Tt4 + (2.24e-8)*Tt4^2 + (C^2)*(exp(C)/(exp(C)-1)^2));
Cp4_fuel = R*(4.47659 + 8.01994e-3*Tt4 - 1.873e-6*(Tt4^2));

% alpha = (h4_air - h3)/h4_fuel;
alpha = Cp4_air*(Tt4-Tt3)/h4_fuel;
%alpha = 0.5*((Cp4_fuel-1)+sqrt((1-Cp4_fuel)^2 + 4*Cp4_air*((Tt4 - Tt3)/h4_fuel)));
alpha_ = alpha*(1-x); %effective richness

h4 = (h4_air+alpha*h4_fuel)/(R+R*alpha);

Tt = [Tt; theta_4];
Pt = [Pt; delta_4];
Ht = [Ht; h4/h0];
% Bleeding

h5_bleed = h4 - (h3 - h2)/((1 - alpha)*(1 - x));

% Turbine:

h5 = (h5_bleed*(1 + alpha)*(1 - x) +h3*x)/((1 + alpha)*(1 - x) + x);
TauT = 1 - (theta_0/theta_4)*(TauC - 1);
theta_5 = TauT*theta_4; % Tt5/T0
pi_t = TauT^(gamma/(eta_tp*(gamma - 1)));
delta_5 = pi_t*delta_4; % Pt5/P0

Tt = [Tt; theta_5];
Pt = [Pt; delta_5];
Ht = [Ht; h5/h0];
% Nozzle:

theta_9 = theta_5; % Tt9/T0
Tt9 = theta_9*T0;
delta_9 = (1 - epsilon_i)*(1 - epsilon_b)*(1 - epsilon_n)*theta_0^(gamma/(gamma - 1))*TauC^(gamma*eta_cp/(gamma - 1))*(1 - (theta_0/theta_4)*(TauC - 1))^(gamma/(eta_tp*(gamma - 1))); % Pt9/P0
Pt9 = delta_9*P0;
B = 3090/Tt9;
h9 = R*(3.5*Tt9 - (Tt9^2)*1.4e-5 + (Tt9^3)*7.467e-9 + 3090/(B-1));
epsilon_T = ((gamma - 1)/gamma)*(epsilon_i + epsilon_b + epsilon_n);
v9 = real(fi*sqrt(theta_4*(1-(1+epsilon_T)/((theta_0*TauC^eta_cp)*(1 - (theta_4/theta_0)*(TauC - 1))^((1 - eta_tp)/eta_tp))) - theta_0*(TauC - 1)));
%V9 = sqrt(2*Cp_air*T0)*sqrt(theta_4 - (TauC*theta_0 - theta_0) - 1);
V9 = sqrt(2*Cp_air*T0)*v9;

Tt = [Tt; theta_9];
Pt = [Pt; delta_9];
Ht = [Ht; h9/h0];
% Performances:

S_T = V9*(1+alpha_) - V0; % Specific thrust
Eta_Overall = (S_T*V0)/(alpha_*h4_fuel);  % Overall efficiency
C_TS = alpha_/S_T; % Thrust specific fuel consumption

end