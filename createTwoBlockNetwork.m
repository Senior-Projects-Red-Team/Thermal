function [network] = createTwoBlockNetwork(constants)

Greenhouse_Temp = 25 + 273.15;
Regolith_Temp = 125;
h_conv = 25; % Estimate from research may need to be fixed

%% Heater Properties
% Water Heater 1 -> 1
% Electronics Heater 2 -> 2
% Life Support Heaters (2x) 3:4 -> 3:4
Rs_Heaters = zeros(4,1);
Cps_Heaters = zeros(4,1) + constants.greenhouse.Heater_cp;
Masses_Heaters = zeros(4,1) + constants.greenhouse.HeaterMass;
Qs_Heaters = Masses_Heaters.*Cps_Heaters.*Greenhouse_Temp;

Rs_Heaters(1) = 0.5; % Estimate fix later.
Rs_Heaters(2:4) = Rs_Heaters(2:4) + 1/(h_conv*constants.greenhouse.HeaterArea);

%% Structure Properties
Rs_Structure = zeros(4,1);
Qs_Structure = zeros(4,1);
Cps_Structure = zeros(4,1);
Masses_Structure = zeros(4,1);

Masses_Structure(1) = 3.1415*constants.greenhouse.radius_inner^2*constants.greenhouse.plate_thickness * constants.greenhouse.structure_density;
Masses_Structure(2) = constants.greenhouse.mass_inner;
Masses_Structure(3) = 0.01; % Estimate for testing, will update later
Masses_Structure(4) = constants.greenhouse.mass_outer;

Cps_Structure(1:2) = constants.greenhouse.structure_cp;
Cps_Structure(3) = 1.8; % kj/kg Cp of PLA estimate
Cps_Structure(4) = constants.greenhouse.structure_cp;

Qs_Structure = Masses_Structure.*Cps_Structure.*Greenhouse_Temp;

% Rs_Structure(1) = constants.greenhouse.plate_thickness / (3.1415*constants.greenhouse.radius_inner^2 * constants.greenhouse.structure_k);
Rs_Structure(1) = 0.001; % High value for ode45
Rs_Structure(2) = 0.328; % Hand Calculated
Rs_Structure(3) = 4.196; % Hand Calculated
Rs_Structure(4) = 0.001; % High calue for ode45

% Plate = 1; -> 5
% innerWall = 2; -> 6
% Bridges = 3; -> 7
% outerWall = 4; -> 8
% Conv to Plate = 9 % No mass props
% Conv to Walls = 10 % No mass props

%% Internal Properties;
% Water = 1 (Note here R is from the water to the Life Support Module) -> 11
% Electronics Block = 2 -> 12
% Life Support Block = 3 -> 13

Masses_Internals = zeros(3,1);
Cps_Internals = zeros(3,1);
Rs_Internals = zeros(3,1);

Masses_Internals(1) = constants.greenhouse.water_mass;
Masses_Internals(2:3) = Masses_Internals(2:3) + (constants.greenhouse.VolWorking/2) * constants.greenhouse.atm_density; % Assumes each module has 1/2 of air.

Cps_Internals(1) = constants.greenhouse.water_cp;
Cps_Internals(2:3) = Cps_Internals(2:3) + constants.greenhouse.atm_cp;

Qs_Internals = Masses_Internals.*Cps_Internals.*Greenhouse_Temp;

Rs_Internals(1) = 0.5; % Estimate of R from water tank to life support

% Note that Rs_Internals(2:3) are 0, as they are unused.

%% Convective Rs

R_convWall = 1/(h_conv.*(constants.greenhouse.AreaWorking/2)); % Working area is halved as each of top and bottome represent half the area of the wall
R_convPlate = 1/(h_conv.*3.1415*constants.greenhouse.radius_inner^2);

%% Moon Slices
% Now we create concentric slices of lunar regolith
slice_thickness = 0.005;
depth = 2; % Depth we expect heat to penetrate (Determined analytically)
num_slices = depth/slice_thickness;

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
slice_Cps = zeros(length(slice_Rs),1) + constants.regolith.cp/1000;

%% Packing

network.Masses = [Masses_Heaters; Masses_Structure; 1;1; Masses_Internals; slice_masses];
network.Cps = [Cps_Heaters; Cps_Structure; 1;1; Cps_Internals; slice_Cps];
network.Qs = [Qs_Heaters; Qs_Structure; 1;1; Qs_Internals; slice_heats];
network.Rs = [Rs_Heaters; Rs_Structure; R_convPlate;R_convWall; Rs_Internals; slice_Rs];
end