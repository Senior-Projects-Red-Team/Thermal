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

Q_gen = 7;
Q_Lights = 10;
heatersMax = 15;

%[Q_dots, heaters, powerUse] = twoBlockModel(Q_0, Masses, Cps, Rs, Q_gen, Q_Lights, heatersMax, 0);

tic();
[ts,Qs] = ode45(@(t,Qs) twoBlockModel(Qs, Masses, Cps, Rs, Q_gen, Q_Lights, heatersMax, t),[0,15*24*3600],Q_0);
toc()


% Cutting the dataset down a bit
Compression_Factor = 1;
ts_small = ts(1:Compression_Factor:end);
Qs_small = Qs(1:Compression_Factor:end,:);

Qdots = zeros(length(ts_small),length(Rs));
heating = zeros(length(ts_small),4);
powerUse = zeros(length(ts_small),1);

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


% Smoothing

for k = 1:60:length(ts_small)-60
    heating(k:k+59,1) = zeros(60,1) + mean(heating(k:k+59,1));
    heating(k:k+59,2) = zeros(60,1) + mean(heating(k:k+59,2));
    heating(k:k+59,3) = zeros(60,1) + mean(heating(k:k+59,3));
    heating(k:k+59,4) = zeros(60,1) + mean(heating(k:k+59,4));
    powerUse(k:k+59) = zeros(60,1) + mean(powerUse(k:k+59));
end

totalHeating = 0.*ts_small;
for k = 1:length(ts_small)
    totalHeating(k) = sum(heating(k,:));
end

%% Plotting

figure()
hold on
title("Internal Temperatures Over Time")
plot(ts_small/(24*3600), Ts(:,13)-273.15, Color="g")
plot(ts_small/(24*3600), Ts(:,12)-273.15, Color="r")
%plot(ts_small/(24*3600), Ts(:,14)-273.15, Color="b")
ylim([20,30])
yline(22, LineStyle=":")
yline(28, LineStyle=":")
legend(["Life Support","Electronics","Tolerance"])
ylabel("Temperature (^oC)")
xlabel("Time (Days)")
hold off

figure()
hold on
title("Relevant Temperatures Over Time")
plot(ts_small/(24*3600), Ts(:,13)-273.15, Color="g")
plot(ts_small/(24*3600), Ts(:,12)-273.15, Color="r")
plot(ts_small/(24*3600), Ts(:,14)-273.15, Color="b")
yline(22, LineStyle=":")
yline(28, LineStyle=":")
legend(["Life Support","Electronics","Regolith", "Tolerance"])
ylabel("Temperature (^oC)")
xlabel("Time (Days)")
hold off

colors = ["r","g","b","k"];
figure()
hold on
title("Heater Wattages over Time")
for i = 1:4
    plot(ts_small/(24*3600), heating(:,i), Color=colors(i), LineWidth=2)
end
legend(["Water", "Electronics", "Life Support(1)", "Life Support (2)"])
ylabel("Power Consumption (Watts)")
xlabel("Time (Days)")
hold off

figure()
hold on
title("Heater Temps over Time")
for i = 1:4
    plot(ts_small/(24*3600), Ts(:,i)-273, Color=colors(i), LineWidth=2)
end
legend(["Water", "Electronics", "Life Support(1)", "Life Support (2)"])
ylabel("Temperature (^oC)")
xlabel("Time (Days)")
hold off

% figure()
% hold on
% title("Life Support Module Temp over Time")
% plot(ts_small/3600, Ts(:,13)-273.15)
% hold off
% 
% figure()
% hold on
% title("Q_{Dot} of Life Support over Time")
% plot(ts_small/3600, Qdots(:,13))
% hold off
% 
% figure()
% hold on
% title("Electronics Temp over Time")
% plot(ts_small/3600, Ts(:,12)-273.15)
% hold off
% 
% figure()
% hold on
% title("Q_{Dot} of Electronics over Time")
% plot(ts_small/3600, Qdots(:,12))
% hold off
% 
% figure()
% hold on
% title("Intermediate Plate Temp over Time")
% plot(ts_small/3600, Ts(:,5)-273.15)
% hold off
% 
% figure()
% hold on
% title("Q_{Dot} of Plate over Time")
% plot(ts_small/3600, Qdots(:,5))
% hold off
% 
figure()
hold on
title("Power Use over Time")
plot(ts_small/(24*3600), powerUse, Color="b")
plot(ts_small/(24*3600), totalHeating, Color="r")
legend("Total", "Heaters")
ylabel("Power Consumption (Watts)")
xlabel("Time (Days)")
hold off

%% Saving Data
Recompression_Factor = 20;
Qs_smaller = Qs_small(1:Recompression_Factor:end,1:30); % Truncating for saving memory when saving
Temps_smaller = Ts(1:Recompression_Factor:end,1:30); % Truncating for saving memory when saving
times_smaller = ts(1:Recompression_Factor:end);

filename = "twoBlock.mat";
clear Q_dots_rel Qs network Qdots Qs_small Ts
save(filename)