clc;
clear;
close all;

%% Running the Two Block Model

constants = loadConstants();
network = createTwoBlockNetwork(constants);

Masses = network.Masses;
Cps = network.Cps;
Q_0 = network.Qs;
Rs = network.Rs;
T_0 = Q_0./(Masses.*Cps);

Q_gen = 10;
Q_Lights = 10;
heatersMax = 15;

[Q_dots, heaters, powerUse] = twoBlockModel(Q_0, Masses, Cps, Rs, Q_gen,Q_Lights, heatersMax, 0);

