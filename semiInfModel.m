function semiInfModel(Q_dot, constants)

heat_flux = @(k, T_s, T_init, alpha, t) (k.*(T_s-T_init))./((3.1415.*alpha.*t).^0.5);

delta_t = 300; % s
numSteps = (15*24*60*60)/delta_t + 1;
time_vect = transpose(linspace(0,15*24*60*60,numSteps));

temp_trials = 500;

regolith_temp = 120; % Kelvin

Surface_temps = linspace(regolith_temp,constants.greenhouse.max_temp, temp_trials);

R_trials = 100;

% R_range = linspace(0,3,R_trials);
% Temps_from_R = constants.greenhouse.max_temp-(R_range.*Q_dot);

R_vals = (constants.greenhouse.max_temp - Surface_temps)/Q_dot;

heat_worst = heat_flux(constants.regolith.conductivity_alt, Surface_temps, regolith_temp, constants.regolith.diffusivity_alt, time_vect);
heat_best = heat_flux(constants.regolith.conductivity, Surface_temps, regolith_temp, constants.regolith.diffusivity, time_vect);

energy_loss = zeros(1,temp_trials);
for i = 1:temp_trials
    energy_loss(i) = sum(heat_worst(2:end,i))*delta_t;
end

net_heat = Q_dot - heat_worst;

net_energy = zeros(numSteps,temp_trials);
for i = 2:numSteps
    for j = 1:temp_trials
    net_energy(i,j) = sum(net_heat(2:i,j));
    end
end

legend_arr = zeros(1,10);

figure()
hold on
for i = 0:9
    legend_arr(i+1) = Surface_temps((temp_trials/10)*i+1);
    plot(time_vect/(24*3600), heat_worst(:,(temp_trials/10)*i+1)./1000)
end
yline(Q_dot/1000);
legend(num2str(transpose(legend_arr)))
title("Heat Flow through Regolith with specified surface Temps and " + Q_dot+ " (W) heat gen")
ylabel("Heat Flow (kW)")
xlabel("Time (Days)")
hold off

figure()
hold on
for i = 0:9
    legend_arr(i+1) = Surface_temps((temp_trials/10)*i+1);
    plot(time_vect/(24*3600), net_energy(:,(temp_trials/10)*i+1)./1000)
end
legend(num2str(transpose(legend_arr)))
yline(0);
title("Net Greenhouse Energy with specified surface Temps and " + Q_dot+ " (W) heat gen")
ylabel("Net Internal Energy (kJ)")
xlabel("Time (Days)")
hold off

figure()
hold on
for i = 0:9
    legend_arr(i+1) = R_vals((temp_trials/10)*i+1);
    plot(time_vect/(24*3600), heat_worst(:,(temp_trials/10)*i+1)./1000)
end
yline(Q_dot/1000);
legend(num2str(transpose(legend_arr)))
title("Heat Flow through Regolith with specified Heat Resistance and " + Q_dot+ " (W) heat gen")
ylabel("Heat Flow (kW)")
xlabel("Time (Days)")
hold off

figure()
hold on
for i = 0:9
    legend_arr(i+1) = R_vals((temp_trials/10)*i+1);
    plot(time_vect/(24*3600), (net_energy(:,(temp_trials/10)*i+1) + (constants.greenhouse.HeatMax+constants.greenhouse.HeatMin)*1000/2)./1000)
end
yline(constants.greenhouse.HeatMin,LineStyle="--");
yline(constants.greenhouse.HeatMax,LineStyle="--");
title("Net Greenhouse Energy with specified Heat Resistance and "  + Q_dot+ " (W) heat gen")
legend(num2str(transpose(legend_arr)))
ylabel("Total Internal Energy (kJ)")
xlabel("Time (Days)")
hold off

%% Finding target values using this

% R_targ = 4; % Inspection of the graph reveals that this will allow 0 heaters.
% 
% ks = linspace(0.001, 0.1, 1001); % A Range of somewhat reasonable k values
% 
% Ls = R_targ*ks*constants.greenhouse.AreaWorking;
% 
% figure()
% hold on
% title("Wall thickness with Thermal Conductivities to hit R_{targ}")
% plot(ks, Ls)
% xline(0.005)
% xline(0.02)
% xline(0.04)
% legend("Thickness vs k", "Vacuum Insulation Panels", "Urethane Foam", "Mineral Wools")
% xlabel("Thermal Conductivity (W/mK)")
% ylabel("")
% hold off

end