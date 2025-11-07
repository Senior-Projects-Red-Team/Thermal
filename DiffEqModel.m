clc;
clear;
close all;

%% Background.

% This is currently the main driver file.

%% Setup
tic();

% This loads in a bunch of relevant information about the greenhouse spec's
constants = loadConstants();
network = createLinNetwork(constants);


%% Running Model

q_gen_base = 10;
q_gen_lights = 10;
powerMax = 40;
heatPumpMax = 0;
heatersOn = 1;

% [ts_actual, Qs_actual, Q_dots, heating, powerUse] = oneDHeatControlledModel(q_gen_base, q_gen_lights, powerMax, heatPumpMax, heatersOn, network);
% 
% powerMax = 50;
% 
% [ts_actual_1, Qs_actual_1, Q_dots_1, heating_1, powerUse_1] = oneDHeatControlledModel(q_gen_base, q_gen_lights, powerMax, heatPumpMax, heatersOn, network);

powerMax = 56;

[ts_actual_2, Qs_actual_2, Q_dots_2, heating_2, powerUse_2] = oneDHeatControlledModel(q_gen_base, q_gen_lights, powerMax, heatPumpMax, heatersOn, network);

Temps_actual_2 = zeros(length(ts_actual_2),10);
for j = 1:10
    Temps_actual_2(:,j) = Qs_actual_2(:,j)./(network.masses(j).*network.cps(j))-273.15;
end

figure()
hold on
plot(ts_actual_2,Temps_actual_2(:,1))
hold off

figure()
hold on
% plot(ts_actual./(24*3600),Qs_actual(:,1)./(network.masses(1)*network.cps(1))-273.15, Color="b", LineStyle="--")
% plot(ts_actual_1./(24*3600),Qs_actual_1(:,1)./(network.masses(1)*network.cps(1))-273.15, Color="b", LineStyle=":")
plot(ts_actual_2./(3600),Qs_actual_2(:,2)./(network.masses(2)*network.cps(2))-273.15, Color="b")
plot(ts_actual_2./(3600),Qs_actual_2(:,6)./(network.masses(6)*network.cps(6))-273.15, Color="r")
% 
yline(28)
yline(22)
title("Active heating, 2 cm Vaccuum Gap Insulation (R= 2.17 W/K)")
% legend(["40 W", "50 W ", "56W" , "Lunar Surface (56 W)", "Viable Range"], Location="southeast")
% 
ylabel("Temperature (C)")
xlabel("Mission Time (Days)")
hold off

figure()
hold on
grid on;
plot(ts_actual_2(1:end)./(24*3600), heating_2(1:end), Color="r")
plot(ts_actual_2(1:end)./(24*3600), powerUse_2(1:end), Color="b")
title("Power Usage with " + powerMax + " maximum power draw.")
ylim([0,powerMax*1.1])
ylabel("Power (Watts)")
xlabel("Mission Time (Days)")
legend(["Heaters", "Total Power Use"])
hold off

