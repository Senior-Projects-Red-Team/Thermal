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

Q_gen = 0;
Q_Lights = 0;
heatersMax = 15;

%[Q_dots, heaters, powerUse] = twoBlockModel(Q_0, Masses, Cps, Rs, Q_gen, Q_Lights, heatersMax, 0);

tic();
[ts,Qs] = ode45(@(t,Qs) twoBlockModel(Qs, Masses, Cps, Rs, Q_gen, Q_Lights, heatersMax, t),[0,6*24*3600],Q_0);
toc()



ts_small = ts(1:20:end);
Qs_small = Qs(1:20:end,:);

Qdots = zeros(length(ts_small),113);
heating = zeros(length(ts_small),4);
powerUse = zeros(length(ts_small),1);
Total_Heat_Flow = zeros(length(ts_small),1);

for i = 1:length(ts_small)
    [Qdots(i,:), heating(i,:), powerUse(i)] = twoBlockModel(transpose(Qs_small(i,:)),Masses, Cps, Rs, Q_gen, Q_Lights, heatersMax, ts_small(i));
end

Qdots = Qdots.*1000;




Q_dots_rel = Qdots(:,11:13);

%% Unpacking
Ts = 0.*Qs_small;
for i = 1:20
Ts(:,i) = Qs_small(:,i)./(Masses(i)*Cps(i));
end



figure()
hold on
title("Life Support Module Temp over Time")
plot(ts_small/3600, Ts(:,13)-273.15)
hold off

figure()
hold on
title("Q_{Dot} of Life Support over Time")
plot(ts_small/3600, Qdots(:,13))
hold off

figure()
hold on
title("Electronics Temp over Time")
plot(ts_small/3600, Ts(:,12)-273.15)
hold off

figure()
hold on
title("Q_{Dot} of Electronics over Time")
plot(ts_small/3600, Qdots(:,12))
hold off

figure()
hold on
title("Intermediate Plate Temp over Time")
plot(ts_small/3600, Ts(:,5)-273.15)
hold off

figure()
hold on
title("Q_{Dot} of Plate over Time")
plot(ts_small/3600, Qdots(:,5))
hold off

figure()
hold on
title("Power Use over Time")
plot(ts_small/3600, powerUse)

hold off

