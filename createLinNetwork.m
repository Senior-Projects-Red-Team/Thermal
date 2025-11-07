function [network] = createLinNetwork(constants)

% Model

% Here we're doing some setup for the DiffEq Model. Calculating thermal
% masses, resistances, and initial temperatures of the greenhouse
% components and of slices of lunar regolith.

h_conv_low = 5; % W/m^2*K Based on general range found via brief googling for viable range
h_conv_high = 25;

median_temp = (constants.greenhouse.max_temp + constants.greenhouse.min_temp)/2;

Internal.Heat = constants.greenhouse.VolWorking * constants.greenhouse.atm_density * constants.greenhouse.atm_cp * median_temp;
%Internal.Heat = Internal.Heat + (constants.greenhouse.water_mass*constants.greenhouse.water_cp*median_temp);

Heater_Heat = constants.greenhouse.Heater_cp*constants.greenhouse.HeaterMass * median_temp;

Plate1.mass = 6.27/2; % Rough Estimate, will fix later.
Plate1.Heat = constants.greenhouse.structure_cp*median_temp * Plate1.mass;
Plate1.L = 0.5e-3; % m

Ribs.mass = 0.5; % Rough Estimate, will fix later.
Ribs.Heat = constants.greenhouse.structure_cp*median_temp * Ribs.mass;
Ribs.A = 4e-5; % m^2
Ribs.L = 20e-3; % m

Plate2.mass = 6.27/2; % Rough Estimate, will fix later.
Plate2.Heat = constants.greenhouse.structure_cp*median_temp * Plate2.mass;
Plate2.L = 1.67e-3; % m

Total_Heat = Internal.Heat + Plate1.Heat + Plate2.Heat + Ribs.Heat;

h_conv = h_conv_high;
R_conv = 1/(h_conv*constants.greenhouse.AreaWorking);

R_1 = Plate1.L / (constants.greenhouse.structure_k*constants.greenhouse.AreaWorking); % Rounded to 0.001 for ode45

R_2 = Ribs.L / (constants.greenhouse.structure_k*Ribs.A);

R_3 = Plate2.L / (constants.greenhouse.structure_k*constants.greenhouse.AreaWorking); % Rounded to 0.001 for ode45

R_wall = R_1 + R_2 + R_3;
L_wall = Plate1.L + Plate2.L + Ribs.L;
k_eff_wall = L_wall / (R_wall * constants.greenhouse.AreaWorking);

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

greenhouse_mass = constants.greenhouse.atm_mass;
%greenhouse_mass = greenhouse_mass + constants.greenhouse.water_mass;
network.Heats = [Heater_Heat; Internal.Heat; Plate1.Heat; Ribs.Heat; Plate2.Heat; slice_heats];
network.masses = [constants.greenhouse.HeaterMass; greenhouse_mass; Plate1.mass; Ribs.mass; Plate2.mass; slice_masses];

greenhouse_cp = constants.greenhouse.atm_cp;
%greenhouse_cp = (constants.greenhouse.atm_mass*constants.greenhouse.atm_cp + constants.greenhouse.water_cp*constants.greenhouse.water_mass)/(constants.greenhouse.atm_mass + constants.greenhouse.water_mass);
structure_cps = ones(3,1).*constants.greenhouse.structure_cp;
slice_cps = ones(num_slices,1).*(constants.regolith.cp/1000);

network.cps = [constants.greenhouse.Heater_cp;greenhouse_cp; structure_cps; slice_cps];

network.Rs = [constants.greenhouse.Heater_R;R_conv; 0.001; R_2; 0.001; slice_Rs];

end

