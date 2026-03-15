clc
clear
close all

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

vars = ["TC0","TC1","TC2","TC3","TMP1","TMP2","Pin9","Pin10","Pin11","Pin6"];

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
DIT_len = height(dryIceTest_data);
DIT_times = linspace(0,10*DIT_len, DIT_len+1);


figure()
hold on
title("Initial Dry Ice")
for i = 1:6
    to_plot = [dryIceTest_data{1,vars(i)};dryIceTest_data{:,vars(i)}]+273.15;
    if i <= 4
        to_plot = to_plot.*correction;
    end
    smoothed = smoothdata(to_plot, 'movmean', 5);
    plot(DIT_times, to_plot-273.15)
    legend(vars(1:6), Location="southeast")
    ylim( [-40,35] )
end

CH_len = height(coolingHeater_data);
CH_times = linspace(0,10*CH_len, CH_len+1);

figure()
hold on
title("coolingHeater")
for i = 1:6
    to_plot = [coolingHeater_data{1,vars(i)};coolingHeater_data{:,vars(i)}]+273.15;
    if i <= 4
        to_plot = to_plot.*correction;
    end
    smoothed = smoothdata(to_plot, 'movmean', 5);
    plot(CH_times, to_plot-273.15)
    legend(vars(1:6), Location="southeast")
    ylim( [-40,35] )
end

[DIT_time, DIT_mean_interior_temp, DIT_mean_inner_wall, DIT_mean_outer_wall, DIT_total_power, DIT_derivative, DIT_R_conv, DIT_R_bridge, DIT_R_foam] = plotCalcs(dryIceTest_data, "DIT", correction, 1);
[CH_time, CH_mean_interior_temp, CH_mean_inner_wall, CH_mean_outer_wall, CH_total_power, CH_derivative, CH_R_conv, CH_R_bridge, CH_R_foam] = plotCalcs(coolingHeater_data, "CH", correction, 1);
[RT_time, RT_mean_interior_temp, RT_mean_inner_wall, RT_mean_outer_wall, RT_total_power, RT_derivative, RT_R_conv, RT_R_bridge, RT_R_foam] = plotCalcs(roomTemp_data, "RT", correction, 1);

DIT_Steady_R_conv = DIT_R_conv(end-20:end);
DIT_Steady_R_bridge = DIT_R_bridge(end-20:end);
DIT_Steady_R_foam = DIT_R_foam(end-20:end);

CH_Steady_R_conv = CH_R_conv(end-20:end);
CH_Steady_R_bridge = CH_R_bridge(end-20:end);
CH_Steady_R_foam = CH_R_foam(end-20:end);

mean_R_conv = mean((DIT_Steady_R_conv + CH_Steady_R_conv)./2);
mean_R_bridge = mean((DIT_Steady_R_bridge + CH_Steady_R_bridge)./2);
mean_R_foam = mean((DIT_Steady_R_foam + CH_Steady_R_foam)./2);
