%% Turbojet Simulator

%% Inputs
M0 = 2; % Mach 0
h = 11e3; % Altitude in meters
T0 = 216.5; % Temperature 0 in Kelvin
Tt4 = 1373; % Total Temperature 4 in Kelvin
TauC = 2; % Compressor T ratio Tt3/Tt2
epsilon_i = 0.075; epsilon_b = 0.06; epsilon_n = 0.02; % Pressure loss coefficients
eta_cp = 0.88; eta_tp = 0.93; % Politropic efficiencies
fi = 0.98; % Nozzle velocity coefficient
x = 0.07; % Turnine cooling bleed

%% Simulation
