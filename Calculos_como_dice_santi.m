clear
clear all
clc
%% Declaramos Variables

M0 = 2;
T0 = 216.5;
P0 = 22632;
Tt4 = 1374;
TET = Tt4;
tau_c = 2;
epsilon_i = 0.075;
epsilon_b = 0.06;
epsilon_n = 0.02;
nu_cp = 0.88;
nu_tp = 0.93;
fi_IsentropicVelocity = 0.98;
x = 0.07;
r = 287.15;
gamma = 1.4;


%% 0: Free Stream:

V0 = M0*sqrt(gamma*r*T0);

Tt0 = T0 * (1+ ((gamma-1)/2)*M0^2);
Pt0 = P0 * (1+ ((gamma-1)/2)*M0^2)^(gamma/(gamma-1));


% Tt0 -> Taules:
Prt0 = Interpolar(380,390,Tt0,3.17,3.47);
ht0 = Interpolar(380,390,Tt0,380.58,390.68);

%% 0 -> 2: Intake

ht2 = ht0;
Pt2 = Pt0*(1-epsilon_i);
Tt2 = Tt0;
Prt2 = Prt0;

%% 3: Compressor

Tt3 = Tt2 * tau_c;

%Tt3 -> Taules:
Prt3 = Interpolar(770,780,Tt3,41.16,43.23);
ht3 = Interpolar(770,780,Tt3,788.87,799.79);

Pt3 = Pt2*((Prt3/Prt2)^nu_cp);

%% 4: Burner

Pt4 = Pt3*(1-epsilon_b);

% Tt4 = TET -> Taules:
ht4 = Interpolar(1370,1380,Tt4,1479.16,1491.13);
Prt4 = Interpolar(1370,1380,Tt4,410.3,423);

alfa = (ht4-ht3)/h_fuel(Tt4);

% Tt4 = TET i alfa -> Taules:
A1 = 0.0169;
A2 = 0.0338;
A = alfa;

% Eje X
B1 = 1370;
B2 = 1380;
B = 1374;

C11 = 1512.48;
C12 = 1524.84;
C21 = 1544.42;
C22 = 1557.45;

ht4 = (((B2-B)/(B2-B1))*C11 + ((B-B1)/(B2-B1))*C12)*((A2-A)/(A2-A1)) + (((B2-B)/(B2-B1))*C21 + ((B-B1)/(B2-B1))*C22)*((A-A1)/(A2-A1));

A1 = 0.0169;
A2 = 0.0338;
A = alfa;

% Eje X
B1 = 1370;
B2 = 1380;
B = 1374;

C11 = 471.3;
C12 = 486.3;
C21 = 538.8;
C22 = 556.4;

Prt4 = (((B2-B)/(B2-B1))*C11 + ((B-B1)/(B2-B1))*C12)*((A2-A)/(A2-A1)) + (((B2-B)/(B2-B1))*C21 + ((B-B1)/(B2-B1))*C22)*((A-A1)/(A2-A1));

%% 5: Turbine

ht5 = ht4 - (ht3-ht2)/((1+alfa)*(1-x)); %ht5 del gas sin mezclar con bleed
ht5 = (ht5*(1+alfa)*(1-x) + ht3*x)/((1+alfa)*(1-x)+x); % ht5 mezclando el gas del compresor con el bleed

% ht5 -> Taules:
T1 = Interpolar(1054.22,1065.94,ht5,990,1000); %interpolamos para sacar la Temp total de ht5 para alfa = 0.0169
T2 = Interpolar(1061.41,1073.42,ht5,980,990); %interpolamos para sacar la Temp total de ht5 para alfa = 0.0338
Tt5 = Interpolar(0.0169,0.0338,alfa,T1,T2);

A1 = 0.0169;
A2 = 0.0338;
A = alfa;

% Eje X
B1 = 990;
B2 = 1000;
B = Tt5;

C11 = 120.68;
C12 = 125.74;
C21 = 132.83;
C22 = 138.54;

Prt5 = (((B2-B)/(B2-B1))*C11 + ((B-B1)/(B2-B1))*C12)*((A2-A)/(A2-A1)) + (((B2-B)/(B2-B1))*C21 + ((B-B1)/(B2-B1))*C22)*((A-A1)/(A2-A1));

Pt5 = Pt4*(Prt5/Prt4)^(1/nu_tp);

%% 9: Nozzle

Pt9 = Pt5*(1-epsilon_n);
Tt9 = Tt5;
ht9 = ht5;

A1 = 0.0169;
A2 = 0.0338;
A = alfa;

% Eje X
B1 = 990;
B2 = 1000;
B = Tt9;

C11 = 120.68;
C12 = 125.74;
C21 = 132.83;
C22 = 138.54;

Prt9 = (((B2-B)/(B2-B1))*C11 + ((B-B1)/(B2-B1))*C12)*((A2-A)/(A2-A1)) + (((B2-B)/(B2-B1))*C21 + ((B-B1)/(B2-B1))*C22)*((A-A1)/(A2-A1));

delta0 = Pt9/P0;
pi9 = delta0;
Pr9_is = exp(log(Prt9) - log(pi9));

% Pr9 -> Taules:
a = Interpolar(8.308,8.70,Pr9_is,490,500); %interpolamos para sacar la Temp total de ht5 para alfa = 0
b = Interpolar(8.36,9.02,Pr9_is,490,500); %interpolamos para sacar la Temp total de ht5 para alfa = 0.0159
T9_is = Interpolar(0.0169,0.0338,alfa,a,b);

% T9 -> Taules

A1 = 0.0169;
A2 = 0.0338;
A = alfa;

% Eje X
B1 = 490;
B2 = 500;
B = T9_is;

C11 = 498.69;
C12 = 509.2;
C21 = 504.63;
C22 = 515.35;

h9_is = (((B2-B)/(B2-B1))*C11 + ((B-B1)/(B2-B1))*C12)*((A2-A)/(A2-A1)) + (((B2-B)/(B2-B1))*C21 + ((B-B1)/(B2-B1))*C22)*((A-A1)/(A2-A1));

V9_is = sqrt(2*(ht9-h9_is)*1000);
V9 = fi_IsentropicVelocity*V9_is;

%% Performances:

alfa_ = alfa*(1-x);

SpecificThrust = (1+alfa_)*V9-V0;
OverallEfficiency = (SpecificThrust*V0)/(alfa_*h_fuel(Tt4)*1000);
SpecificFuelConsumption = alfa_/SpecificThrust;






