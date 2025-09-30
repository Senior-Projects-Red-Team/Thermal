clc;
clear;
close all;

%% Background.

% This is currently the main driver file.

%% Setup

% This loads in a bunch of relevant information about the greenhouse spec's
constants = loadConstants();

h_conv_low = 5; % W/m^2*K Based on general range found via brief googling for viable range
h_conv_high = 25;

char_L = constants.greenhouse.VolWorking / constants.greenhouse.AreaWorking;

Biot_Gas_Filled_Cube = @(h, Lc, k) h.*Lc./k;


Biot_greenhouse_high_h = Biot_Gas_Filled_Cube(h_conv_high, char_L, constants.greenhouse.structure_k);
Biot_greenhouse_low_h = Biot_Gas_Filled_Cube(h_conv_low, char_L, constants.greenhouse.structure_k);
% In both cases the Biot number is less than 0.1

tic();

%% Model

% Here we're doing some setup for the DiffEq Model. Calculating thermal
% masses, resistances, and initial temperatures of the greenhouse
% components and of slices of lunar regolith.

median_temp = (constants.greenhouse.max_temp + constants.greenhouse.min_temp)/2;

Internal.Heat = constants.greenhouse.VolWorking * constants.greenhouse.atm_density * constants.greenhouse.atm_cp * median_temp;
Internal.Heat = Internal.Heat + (constants.greenhouse.water_mass*constants.greenhouse.water_cp*median_temp);

Plate1.mass = 5.08/2; % Rough Estimate, will fix later.
Plate1.Heat = constants.greenhouse.structure_cp*median_temp * Plate1.mass;
Plate1.L = 0.9e-3; % m

Ribs.mass = 1.67; % Rough Estimate, will fix later.
Ribs.Heat = constants.greenhouse.structure_cp*median_temp * Ribs.mass;
Ribs.A = 4e-5; % m^2
Ribs.L = 20e-3; % m

Plate2.mass = 5.08/2; % Rough Estimate, will fix later.
Plate2.Heat = constants.greenhouse.structure_cp*median_temp * Plate2.mass;
Plate2.L = 1.3e-3; % m

Total_Heat = Internal.Heat + Plate1.Heat + Plate2.Heat + Ribs.Heat;

h_conv = h_conv_high;
R_conv = 1/(h_conv*constants.greenhouse.AreaWorking);

R_1 = Plate1.L / (constants.greenhouse.structure_k*constants.greenhouse.AreaWorking); % Rounded to 0.001 for ode45

R_2 = Ribs.L / (constants.greenhouse.structure_k*Ribs.A);

R_3 = Plate2.L / (constants.greenhouse.structure_k*constants.greenhouse.AreaWorking); % Rounded to 0.001 for ode45

R_wall = R_1 + R_2 + R_3;
L_wall = Plate1.L + Plate2.L + Ribs.L;
k_eff_wall = L_wall / (R_wall * constants.greenhouse.AreaWorking);

Biot_greenhouse_effective_wall_low_h = Biot_Gas_Filled_Cube(h_conv_low, char_L, k_eff_wall);
Biot_greenhouse_effective_wall_high_h = Biot_Gas_Filled_Cube(h_conv_high, char_L, k_eff_wall);
% Now both of these are much higher than 0.1, which is why we need to
% analyze them seperately, not as a lumped sum.


h_c_Al = 42000; % (Chosen for lowest estimate, and thus worst insulation)
R_contact = 1/(constants.greenhouse.AreaWorking*h_c_Al); % Neglected further on, as near 0

R_total = R_conv + R_1 + R_2 + R_3 + R_contact;

% Now we create concentric slices of lunar regolith
num_slices = 400;
depth = 2; % Depth we expect heat to penetrate (Determined analytically)
slice_thickness = depth / num_slices;

slice_max_area = 4*3.1415927*(((.35/2)+2)^2); % Note that this is spherical
slice_min_area = constants.greenhouse.AreaWorking;

slice_areas = transpose(linspace(slice_min_area,slice_max_area,num_slices));

% slice_areas = zeros(num_slices,1);
% for i = 0:num_slices-1
%     slice_areas(i+1) = 6*((0.35 + 2*(slice_thickness*i))^2);
% end

% slice_areas = ones(num_slices,1)*constants.greenhouse.AreaWorking;


slice_volumes = slice_thickness.*slice_areas;
slice_masses = slice_volumes.*constants.regolith.density;

slice_heats = slice_masses.*(constants.regolith.cp/1000).*constants.regolith.Maxtemp;

slice_Rs = slice_thickness./(constants.regolith.conductivity*slice_areas);


%% Combining elements into vectors


Heats = [Internal.Heat; Plate1.Heat; Ribs.Heat; Plate2.Heat; slice_heats];
masses = [constants.greenhouse.water_mass+constants.greenhouse.atm_mass; Plate1.mass; Ribs.mass; Plate2.mass; slice_masses];

greenhouse_cp = (constants.greenhouse.atm_mass*constants.greenhouse.atm_cp + constants.greenhouse.water_cp*constants.greenhouse.water_mass)/(constants.greenhouse.atm_mass + constants.greenhouse.water_mass);
structure_cps = ones(3,1).*constants.greenhouse.structure_cp;
slice_cps = ones(num_slices,1).*(constants.regolith.cp/1000);

cps = [greenhouse_cp; structure_cps; slice_cps];

Rs = [R_conv; 0.001; R_2; 0.001; slice_Rs];

Temps = Heats./(masses.*cps);

Q_gen = 50;


%[ts,Qs] = ode45(@(t,Qs) oneDWithHeatPump(Qs, masses, cps, Rs, 10, 10, 0, 64, 0, t),[0,15*24*3600],Heats);
q_gen_base = 10;
q_gen_lights = 10;
powerMax = 70;
heatPumpMax = 0;
heatersOn = 1;


% options = odeset('OutputFcn', @(t,y,flag) myOutputFcn(t,y,flag, m0, Cload));


%options = odeset('Events',heatersOnEvent);
[ts,Qs] = ode45(@(t,Qs) oneDWithHeatPump(Qs, masses, cps, Rs, q_gen_base, q_gen_lights, heatPumpMax, powerMax, heatersOn, t),[0,15*24*3600],Heats);
toc()

%Extract extra variables by re-evaluating the ODE function
Q_dots = 0.*Qs;
heating = zeros(length(ts),1);
powerUse = zeros(length(ts),1);
for i = 1:length(ts)
    [Q_dots(i,:), heating(i), powerUse(i)] = oneDWithHeatPump(transpose(Qs(i,:)), masses, cps, Rs, q_gen_base, q_gen_lights, heatPumpMax, powerMax, heatersOn, ts(i));
end
toc()

figure()
hold on
plot(ts./(24*3600),Qs(:,1)./(masses(1)*cps(1))-273.15, Color="b")
plot(ts./(24*3600),Qs(:,5)./(masses(5)*cps(5))-273.15, Color="r")
yline(28)
yline(22)
hold off

figure()
hold on
grid on;
plot(ts(1:1000:end)./(24*3600), heating(1:1000:end))
plot(ts(1:1000:end)./(24*3600), powerUse(1:1000:end))
legend(["heating", "powerUse"])
hold off

total_Power_Use = (sum(powerUse)/15*24*3600)/1000;

% runModelTemps(masses, cps, Rs, Heats, 50,0,0,0,1,0,1);
% toc()
% runModelTemps(masses, cps, Rs, Heats, 50,50,1,0,1,0,2);
% runModelTemps(masses, cps, Rs, Heats, 64,0,0,0,1,0,3);
% runModelTemps(masses, cps, Rs, Heats, 64,64,1,0,1,0,4);
% 
% 
% runModelTemps(masses, cps, Rs, Heats, 70,0,0,0,1,0,5);
% runModelTemps(masses, cps, Rs, Heats, 70,70,1,0,1,0,6);


% runModelHeats(masses, cps, Rs, Heats, 50, 2, constants)

% semiInfModel(50, constants)

toc()

function runModelTemps(Masses, Cps, Rs, Heats, Q_gen, heatPump, heatersBool, wallBool, surfaceBool, slicesBool, figNum)

figure(figNum)
hold on
for i = 1:(15*24)

    if(i == 1)
        [ts,Qs] = ode45(@(t,Qs) oneDWithHeatPump(Qs, Masses, Cps, Rs, Q_gen, heatPump, heatersBool, t),[0,i*3600],Heats);
    else
        [ts,Qs] = ode45(@(t,Qs) oneDWithHeatPump(Qs, Masses, Cps, Rs, Q_gen, heatPump, heatersBool, t),[(i-1)*3600,i*3600],Qs(end,:));
    end
    plot(ts./(24*3600),Qs(:,1)./(Masses(1)*Cps(1))-273.15,Color="b")

    if(wallBool)
        for j = 2:4
            plot(ts./(24*3600),Qs(:,j)./(Masses(j)*Cps(j))-273.15,Color="r",LineStyle="-.")
        end
    end

    if(surfaceBool)
        plot(ts./(24*3600),Qs(:,5)./(Masses(5)*Cps(5))-273.15,Color="r",LineStyle="-.")
    end
    if(slicesBool)
        for j = 5:length(Rs)
            if(mod(j,5)==0)
                plot(ts./(24*3600),Qs(:,j)./(Masses(j)*Cps(j))-273.15,LineStyle=":",Color="g")
            end
        end
    end

    if(i == 1)
        yline(28)
        yline(22)
        title("Temp of Internal Atmosphere with " + Q_gen + " (W) Internal Heat Gen and" + heatPump + "(W) Heat Pump")
        ylabel("C")
        xlabel("Time (Days)")
    end

end
legend(["Internal", "Lunar Surface", "Viable Range"],Location="southeast")
hold off

end

function runModelHeats(Masses, Cps, Rs, Heats, Q_gen, figNum, constants)

fig = figNum+10;

figure(fig)
hold on
for i = 1:(15*24)

    if(i == 1)
        [ts,Qs] = ode45(@(t,Qs) oneDHeatFlowModel(Qs, Masses, Cps, Rs, Q_gen, t),[0,i*3600],Heats);
    else
        [ts,Qs] = ode45(@(t,Qs) oneDHeatFlowModel(Qs, Masses, Cps, Rs, Q_gen, t),[(i-1)*3600,i*3600],Qs(end,:));
    end

    greenhouse_heat = Qs(:,1) + Qs(:,2) + Qs(:,3) + Qs(:,4); 

    plot(ts./(24*3600),greenhouse_heat,Color="b")
end

yline(constants.greenhouse.HeatMin)
yline(constants.greenhouse.HeatMax)
title("Total Greenhouse Heat over time with " + Q_gen + " (W) Internal Heat Gen")
ylabel("Heat (kJ)")
xlabel("Time (Days)")

hold off

end