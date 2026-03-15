clc;
clear;
close all;

%% Parsing data and calculating derived values.


coolingHeater_filename = "CoolTerm Capture (Cooling_Heater_dry_Ice_Test.CoolTermSettings) 2026-02-20 10-55-56-743.txt";
Dry_Ice_Test_filename = "CoolTerm Capture (Dry_Ice_Test_Data.CoolTermSettings) 2026-02-20 08-50-06-140.txt";
roomTemp_filename = "CoolTerm Capture (Room_temp_test.CoolTermSettings) 2026-02-20 12-22-10-833.txt";

coolingHeater_data = Parse_Data(coolingHeater_filename);
dryIceTest_data = Parse_Data(Dry_Ice_Test_filename);
roomTemp_data = Parse_Data(roomTemp_filename);

end_state = roomTemp_data{end,:} + 273.15;
true_end_temp = mean(end_state(5:6));
thermocouple_end_temp = mean(end_state(1:4));
correction = true_end_temp / thermocouple_end_temp;
% 
% RT_len = height(roomTemp_data);
% RT_times = linspace(0,10*RT_len, RT_len+1);

% vars = ["TC0","TC1","TC2","TC3","TMP1","TMP2","Pin9","Pin10","Pin11","Pin6"];

% figure()
% hold on
% title("Reaching Room Temp")
% for i = 1:6
%     to_plot = [roomTemp_data{1,vars(i)};roomTemp_data{:,vars(i)}] + 273.15;
%     if i <= 4
%         to_plot = to_plot.*correction;
%     end
%     smoothed = smoothdata(to_plot, 'movmean', 3);
%     plot(RT_times, smoothed-273.15)
%     legend(vars(1:6), Location="southeast")
%     ylim( [-40,35] )
% end
% 
% DIT_len = height(dryIceTest_data);
% DIT_times = linspace(0,10*DIT_len, DIT_len+1);

% 
% figure()
% hold on
% title("Initial Dry Ice")
% for i = 1:6
%     to_plot = [dryIceTest_data{1,vars(i)};dryIceTest_data{:,vars(i)}]+273.15;
%     if i <= 4
%         to_plot = to_plot.*correction;
%     end
%     smoothed = smoothdata(to_plot, 'movmean', 5);
%     plot(DIT_times, smoothed-273.15)
%     legend(vars(1:6), Location="southeast")
%     ylim( [-40,35] )
% end

% CH_len = height(coolingHeater_data);
% CH_times = linspace(0,10*CH_len, CH_len+1);
% 
% figure()
% hold on
% title("coolingHeater")
% for i = 1:6
%     to_plot = [coolingHeater_data{1,vars(i)};coolingHeater_data{:,vars(i)}]+273.15;
%     if i <= 4
%         to_plot = to_plot.*correction;
%     end
%     smoothed = smoothdata(to_plot, 'movmean', 5);
%     plot(CH_times, smoothed-273.15)
%     legend(vars(1:6), Location="southeast")
%     ylim( [-40,35] )
% end

[DIT_time, DIT_mean_interior_temp, DIT_mean_inner_wall, DIT_mean_outer_wall, DIT_total_power, DIT_derivative, DIT_R_conv, DIT_R_bridge, DIT_R_foam] = plotCalcs(dryIceTest_data, "DIT", correction, 0);
[CH_time, CH_mean_interior_temp, CH_mean_inner_wall, CH_mean_outer_wall, CH_total_power, CH_derivative, CH_R_conv, CH_R_bridge, CH_R_foam] = plotCalcs(coolingHeater_data, "CH", correction, 0);
% plotCalcs(roomTemp_data, "RT", correction)

DIT_Steady_R_conv = DIT_R_conv(end-5:end);
DIT_Steady_R_bridge = DIT_R_bridge(end-5:end);
DIT_Steady_R_foam = DIT_R_foam(end-5:end);

CH_Steady_R_conv = CH_R_conv(end-5:end);
CH_Steady_R_bridge = CH_R_bridge(end-5:end);
CH_Steady_R_foam = CH_R_foam(end-5:end);

% mean_R_conv = mean((DIT_Steady_R_conv + CH_Steady_R_conv)./2);
% mean_R_bridge = mean((DIT_Steady_R_bridge + CH_Steady_R_bridge)./2);
% mean_R_foam = mean((DIT_Steady_R_foam + CH_Steady_R_foam)./2);


%% Running the Two Block Model

constants = loadConstants();
network = createTwoBlockNetwork(constants);

Masses = network.Masses(1:13);
Cps = network.Cps(1:13);
Q_0 = network.Qs(1:13);
Rs = network.Rs(1:13);
T_0 = Q_0./(Masses.*Cps);


Rs(6) = 0.1; % This one is lower, because without a thermal gap the resistance along the inner wall is negligible compared to convection across the gap.

k_air_gap = 0.026;
R_air_gap_cap = 0.01 / (k_air_gap * 3.1415 * 0.17^2);
R_air_gap_wall = log(0.17/0.16)/(2*3.1415*0.36 * k_air_gap);
R_air_gap = (2/R_air_gap_cap + R_air_gap_wall)^(-1);

Rs(7) = R_air_gap; % Testing a lower insulation (Lower than 1.6 results in plants getting too cold, although repositioning water heater may allow lower Rs to be used)
              % Raising the target temp to 27 C allows system to handle
              % insulation as low as 1.3. Implementing Integral  control
              % may also allow for lower Rs(7). Note that lower insulation
              % requires more power.
% Here we tweak the resistance values based on experimental data.
Rs(7) = mean(CH_Steady_R_bridge);
Rs(8) = mean(CH_Steady_R_foam); % This is acting to capture the effect of the foam on the exterior.
Rs(10) = 2.*mean(CH_Steady_R_conv);
%
%Rs(10) = 0.4;
Q_gen = 16*0.5;
Q_Lights = 16*0.67;
heatersMax = 16;
T_ext = -78.5 +273.15;

%[Q_dots, heaters, powerUse] = twoBlockModel(Q_0, Masses, Cps, Rs, Q_gen, Q_Lights, heatersMax, 0);

tic();
[ts,Qs] = ode45(@(t,Qs) twoBlockTestModel(Qs, Masses, Cps, Rs, Q_gen, Q_Lights, heatersMax, T_ext , t),[0,10000],Q_0);
toc()

% Cutting the dataset down a bit
Compression_Factor = 1;
ts_small = ts(1:Compression_Factor:end);
delta_ts_small = max(ts_small)/length(ts_small);
Qs_small = Qs(1:Compression_Factor:end,:);

Qdots = zeros(length(ts_small),length(Rs));
heating = zeros(length(ts_small),4);
powerUse = zeros(length(ts_small),1);

for i = 1:length(ts_small)
    [Qdots(i,:), heating(i,:), powerUse(i)] = twoBlockTestModel(transpose(Qs_small(i,:)),Masses, Cps, Rs, Q_gen, Q_Lights, heatersMax, T_ext , ts_small(i));
end


Qdots = Qdots.*1000;


Q_dots_rel = Qdots(:,11:13);

%% Unpacking
Ts = 0.*Qs_small;
for i = 1:13
Ts(:,i) = Qs_small(:,i)./(Masses(i)*Cps(i));
end


% Smoothing

% for k = 1:60:length(ts_small)-60
%     heating(k:k+59,1) = zeros(60,1) + mean(heating(k:k+59,1));
%     heating(k:k+59,2) = zeros(60,1) + mean(heating(k:k+59,2));
%     heating(k:k+59,3) = zeros(60,1) + mean(heating(k:k+59,3));
%     heating(k:k+59,4) = zeros(60,1) + mean(heating(k:k+59,4));
%     powerUse(k:k+59) = zeros(60,1) + mean(powerUse(k:k+59));
% end
heating(:,1) = smoothdata(heating(:,1), 'movmean', 60);
heating(:,2) = smoothdata(heating(:,2), 'movmean', 60);
heating(:,3) = smoothdata(heating(:,3), 'movmean', 60);
heating(:,4) = smoothdata(heating(:,4), 'movmean', 60);
powerUse = smoothdata(powerUse, 'movmean', 60);


totalHeating = 0.*ts_small;
for k = 1:length(ts_small)
    totalHeating(k) = sum(heating(k,:));
end

%% Plotting

figure()
hold on
title("Internal Temperatures Over Time")
plot(ts_small, Ts(:,13)-273.15, Color="g")
plot(ts_small, Ts(:,12)-273.15, Color="r")
%plot(ts_small/(24*3600), Ts(:,14)-273.15, Color="b")
ylim([0,50])
yline(22, LineStyle=":")
yline(28, LineStyle=":")
legend(["Life Support","Electronics","Tolerance"])
ylabel("Temperature (^oC)")
xlabel("Time (Days)")
hold off

figure()
hold on
title("Relevant Temperatures Over Time")
plot(ts_small, Ts(:,13)-273.15, Color="g")
plot(ts_small, Ts(:,12)-273.15, Color="r")
yline(-80, Color="b")
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
    plot(ts_small, heating(:,i), Color=colors(i), LineWidth=2)
end
legend(["Water", "Electronics", "Life Support(1)", "Life Support (2)"])
ylabel("Power Consumption (Watts)")
xlabel("Time (Days)")
hold off

figure()
hold on
title("Heater Temps over Time")
for i = 1:4
    plot(ts_small, Ts(:,i)-273, Color=colors(i), LineWidth=2)
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
DIT_total_power(1:256) = DIT_total_power(1:256).*0;
figure()
hold on
title("Power Use over Time")
plot(ts_small, powerUse, Color="b")
plot(ts_small, totalHeating, Color="r")
plot(DIT_time+0, DIT_total_power, Color="g")
legend("Total", "Heaters", "Total Experimental")
ylabel("Power Consumption (Watts)")
xlabel("Time (sec)")
hold off

figure()
hold on
title("Internal Temperatures Over Time")
plot(ts_small, Ts(:,13), Color="g")
plot(ts_small, Ts(:,12), Color="r")
plot(ts_small, Ts(:,6))
%plot(ts_small, Ts(:,7))
plot(ts_small, Ts(:,8))
plot(DIT_time, DIT_mean_interior_temp+273.15, LineStyle="--")
plot(DIT_time, DIT_mean_inner_wall+273.15, LineStyle="--")
plot(DIT_time, DIT_mean_outer_wall+273.15, LineStyle="--")
legend(["Life Support","Electronics","Inner Wall", "Outer Wall"])
ylim([0,350])
ylabel("Temperature (^oC)")
xlabel("Time (Days)")
hold off

%% Saving Data
Recompression_Factor = 20;
Qs_smaller = Qs_small(1:Recompression_Factor:end,1:13); % Truncating for saving memory when saving
Temps_smaller = Ts(1:Recompression_Factor:end,1:13); % Truncating for saving memory when saving
times_smaller = ts(1:Recompression_Factor:end);

filename = "twoBlockTest.mat";
clear Q_dots_rel Qs network Qdots Qs_small Ts
save(filename)