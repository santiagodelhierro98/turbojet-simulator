%% Turbojet Simulator
clear all;
close all;
clc;
clear;

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
x = 0.07; % Turbine cooling bleed
R = 287.15; % Gas cte in J/kgK
hf0 = 4.3095e7; % Fuel lower heating value J/kg

%% Simulation
[Tt, Pt, Ht, S_T, Eta_Overall, C_TS, Eta_Thermal, Eta_Propulsive] = Turbojet_Sim_2(M0, T0, P0, Tt4, TauC, epsilon_i, epsilon_b, epsilon_n, eta_cp, eta_tp, fi, x, R, hf0);
fprintf("Psi = " + S_T + "m/s\n");
fprintf("C_t_s = " + 1000*C_TS + "g/kN*s\n");
fprintf("Eta_o = " + 100*Eta_Overall + "%");
fprintf("Eta_th = " + 100*Eta_Thermal + "%");
fprintf("Eta_pr = " + 100*Eta_Propulsive + "%");

figure()
subplot(3, 1, 1)
plot(Tt, '-o')
xlabel('Stage');
ylabel('Total Temperature [K]');
grid on

subplot(3, 1, 2)
plot(Pt, '-o')
xlabel('Stage');
ylabel('Total Pressure [Pa]');
grid on

subplot(3, 1, 3)
plot(Ht, '-o')
xlabel('Stage');
ylabel('Enthalpy [J/kg]');
grid on

S_Tv = [];
Eta_Overallv = [];
C_TSv = [];
Tt4v = [];

n = 1;
M0 = 0;
while (M0 <= 2.5)
    Tt4 = 1000;
    i = 1;
    while (Tt4 <= 1900)
        [Tt, Pt, Ht, S_T, Eta_Overall, C_TS] = Turbojet_Sim_2(M0, T0, P0, Tt4, TauC, epsilon_i, epsilon_b, epsilon_n, eta_cp, eta_tp, fi, x, R, hf0);
        S_Tv(n, i) = S_T;
        if (Eta_Overall < 0)
            Eta_Overall = 0;
        end
        Eta_Overallv(n, i) = Eta_Overall;
        C_TSv(n, i) = C_TS;
        Tt4v(i) = Tt4;
        
        i= i +1;
        Tt4 = Tt4 + 10;
    end
    n = n + 1;
    M0 = M0 + 0.5;
end

figure()
subplot(1,3,1)
plot(Tt4v, S_Tv());
title('Specific Thrust');
xlabel('T_t_4 [K]');
ylabel('Psi [m/s]');
legend('M_0 = 0', 'M_0 = 0.5', 'M_0 = 1', 'M_0 = 1.5', 'M_0 = 2','M_0 = 2.5');
grid on;

subplot(1,3,2)
plot(Tt4v, 100*Eta_Overallv)
title('Overall Efficiency');
xlabel('T_t_4 [K]');
ylabel('Eta_o [%]');
legend('M_0 = 0', 'M_0 = 0.5', 'M_0 = 1', 'M_0 = 1.5', 'M_0 = 2','M_0 = 2.5');
grid on;

subplot(1,3,3)
plot(Tt4v, 1000*C_TSv)
title('Thrust Specific Fuel Consumption');
xlabel('T_t_4 [K]');
ylabel('C_t_s [g/kN*s]');
legend('M_0 = 0', 'M_0 = 0.5', 'M_0 = 1', 'M_0 = 1.5', 'M_0 = 2','M_0 = 2.5');
grid on;