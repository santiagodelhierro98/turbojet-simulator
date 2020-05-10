%% Turbojet Simulator
clear all;
close all;

%% Inputs
% M0 = 2; % Mach 0
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
S_Tv = [];
Eta_Overallv = [];
C_TSv = [];
Tt4v = [];
M0v = [];

Tt4 = 1000;
n = 1;
while(Tt4 <= 1800)
    M0 = 0;
    i = 1;
    while (M0 <= 2.45)
        [S_T, Eta_Overall, C_TS] = Turbojet_Sim(M0, T0, P0, Tt4, TauC, epsilon_i, epsilon_b, epsilon_n, eta_cp, eta_tp, fi, x, R, hf0);
        S_Tv(n, i) = S_T;
        if (Eta_Overall < 0)
            Eta_Overall = 0;
        end
        Eta_Overallv(n, i) = Eta_Overall;
        C_TSv(n, i) = C_TS;
        M0v(i) = M0;
        
        i= i +1;
        M0 = M0 + 0.01;
    end
    
    n = n + 1;
    Tt4 = Tt4 + 400;    
end

subplot(1,3,1)
plot(M0v, S_Tv);
title('Specific Thrust');
xlabel('Mach');
ylabel('Psi [m/s]');
legend('T_t_4 = 1000 K', 'T_t_4 = 1400 K', 'T_t_4 = 1800 K');
grid on;

subplot(1,3,2)
plot(M0v, 100 * Eta_Overallv)
title('Overall Efficiency');
xlabel('Mach');
ylabel('Eta_o [%]');
legend('T_t_4 = 1000 K', 'T_t_4 = 1400 K', 'T_t_4 = 1800 K');
grid on;

subplot(1,3,3)
plot(M0v, 1000 * C_TSv)
title('Thrust Specific Fuel Consumption');
xlabel('Mach');
ylabel('C_t_s [g/kN*s]');
legend('T_t_4 = 1000 K', 'T_t_4 = 1400 K', 'T_t_4 = 1800 K');
grid on;
