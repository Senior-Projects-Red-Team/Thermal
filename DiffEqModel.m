clc;
clear;
close all;

%% Setup

constants = loadConstants();

tic();

%% Model

median_temp = (constants.greenhouse.max_temp + constants.greenhouse.min_temp)/2;

Internal.Heat = constants.greenhouse.VolWorking * constants.greenhouse.atm_density * constants.greenhouse.atm_cp * median_temp;
Internal.Heat = Internal.Heat + (constants.greenhouse.water_mass*constants.greenhouse.water_cp*median_temp);

Plate1.mass = 5.08/2; % Rough Estimate, will fix later.
Plate1.Heat = constants.greenhouse.structure_cp*median_temp * Plate1.mass;
Plate1.L = 0.9e-3; % m

Ribs.mass = 1.67; % Rough Estimate, will fix later.
Ribs.Heat = constants.greenhouse.structure_cp*median_temp * Ribs.mass;
Ribs.A = 4e-5; % m^2
Ribs.L = 15e-3; % m

Plate2.mass = 5.08/2; % Rough Estimate, will fix later.
Plate2.Heat = constants.greenhouse.structure_cp*median_temp * Plate2.mass;
Plate2.L = 1.3e-3; % m

Total_Heat = Internal.Heat + Plate1.Heat + Plate2.Heat + Ribs.Heat;


R_conv = 0.5;

R_1 = Plate1.L / (constants.greenhouse.structure_k*constants.greenhouse.AreaWorking);

R_2 = Ribs.L / (constants.greenhouse.structure_k*Ribs.A);

R_3 = Plate2.L / (constants.greenhouse.structure_k*constants.greenhouse.AreaWorking);


h_c_Al = 42000; % (Chosen for lowest estimate, and thus worst insulation)
R_contact = 1/(constants.greenhouse.AreaWorking*h_c_Al);

R_total = R_conv + R_1 + R_2 + R_3 + R_contact;

% Now we create concentric slices of lunar regolith
num_slices = 200;
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

% runmodel(masses, cps, Rs, Heats, 32,0,1,0,1);
runModelTemps(masses, cps, Rs, Heats, 50,0,1,0,2);
% runmodel(masses, cps, Rs, Heats, 64,0,1,0,3);
% runmodel(masses, cps, Rs, Heats, 96,0,1,0,4);
% runmodel(masses, cps, Rs, Heats, 128,0,1,0,5);

runModelHeats(masses, cps, Rs, Heats, 50, 2, constants)

semiInfModel(50, constants)

toc()

function runModelTemps(Masses, Cps, Rs, Heats, Q_gen, wallBool, surfaceBool, slicesBool, figNum)

figure(figNum)
hold on
for i = 1:(15*24)

    if(i == 1)
        [ts,Qs] = ode45(@(t,Qs) oneDHeatFlowModel(Qs, Masses, Cps, Rs, Q_gen, t),[0,i*3600],Heats);
    else
        [ts,Qs] = ode45(@(t,Qs) oneDHeatFlowModel(Qs, Masses, Cps, Rs, Q_gen, t),[(i-1)*3600,i*3600],Qs(end,:));
    end
    plot(ts./(24*3600),Qs(:,1)./(Masses(1)*Cps(1))-273.15,Color="b")

    if(wallBool)
        for j = 2:4
            plot(ts./(24*3600),Qs(:,j)./(Masses(j)*Cps(j))-273.15,Color="r",LineStyle="-.")
        end
    end

    if(surfaceBool)
        plot(ts./(24*3600),Qs(:,4)./(Masses(4)*Cps(4))-273.15,Color="r",LineStyle="-.")
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
        title(Q_gen)
    end

end
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

hold off

end