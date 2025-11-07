function constants = loadConstants()

constants.regolith.conductivity = 2.95e-2; % 2.95e-2 Watts/m*K (Worst case scenario, revise as needed) % Maybe 7.4 e-4 W/mK @SP Surface 3.4e-3 W/mK @ 1m ??? % Range 0.5e-3 -> 0.8e-3

constants.regolith.Mintemp = 40; % K
constants.regolith.Maxtemp = 125; % K
constants.regolith.density = 1800; % kg/m^3 @ 1m
constants.regolith.cp = 0.840; % kJ/kgK
constants.regolith.diffusivity = constants.regolith.conductivity / (constants.regolith.density * constants.regolith.cp); % m^2/s; Seems way too low.
constants.regolith.diffusivity_alt = 2e-7; % m^2/s
constants.regolith.conductivity_alt = constants.regolith.diffusivity_alt*constants.regolith.density*constants.regolith.cp;% W/m*K

constants.power.battery_mass = 0;
constants.power.battery_cp = 0.8; % kj/kgK 

constants.greenhouse.max_temp = 28 + 273.15; % K
constants.greenhouse.min_temp = 22 + 273.15; % K
constants.greenhouse.atm_cp = 1.005; % kJ/kgK % Assuming similar to normal air for now
constants.greenhouse.structure_cp = 0.921 ; % 0.921 for Aluminum, 0.75 for Carbon Fiber, 1.5 for Kevlar, HDPE 1.3-2.2, all units as kJ/kgK
constants.greenhouse.structure_k = 167; % W/mK Value for aluminum.
constants.greenhouse.structure_density = 2700; % kg/m^3
constants.greenhouse.water_cp = 4.18; % kJ/kgK
constants.greenhouse.atm_density = 1.204; % kg/m^3 Assumes normal air
constants.greenhouse.water_density = 1000; %kg/m^3
constants.greenhouse.water_mass = 0.075; %
constants.greenhouse.structure_mass = 6.75; % Consult with SE and Structures Teams

constants.greenhouse.HeaterMass = 0.2; % kg
constants.greenhouse.Heater_cp = 0.46; % kJ/kgK
constants.greenhouse.Heater_R = 2.3; % W/K
constants.greenhouse.HeaterArea = 0.05*0.1;

constants.greenhouse.radius_inner = 0.179 ; 
constants.greenhouse.radius_outer = 0.189 ;
constants.greenhouse.radius_avg = (constants.greenhouse.radius_inner+constants.greenhouse.radius_outer)/2;
constants.greenhouse.length = 0.36;

constants.greenhouse.plate_thickness = 0.001; % m

constants.greenhouse.thickness_inner = 0.5e-3; % m
constants.greenhouse.vol_inner = 3.1415 * ((constants.greenhouse.radius_inner+constants.greenhouse.thickness_inner)^2-constants.greenhouse.radius_inner^2) * constants.greenhouse.length;
constants.greenhouse.thickness_outer = 1.67e-3; % m
constants.greenhouse.vol_outer = 3.1415 * ((constants.greenhouse.radius_outer+constants.greenhouse.thickness_outer)^2-constants.greenhouse.radius_outer^2) * constants.greenhouse.length;
constants.greenhouse.mass_inner = constants.greenhouse.structure_density * constants.greenhouse.vol_inner;
constants.greenhouse.mass_outer = constants.greenhouse.structure_density * constants.greenhouse.vol_outer;

constants.greenhouse.thickness_cap = 3.4e-3; % m;
constants.greenhouse.vol_cap = 3.1415*constants.greenhouse.radius_outer^2*constants.greenhouse.thickness_cap;
constants.greenhouse.mass_cap = 2*constants.greenhouse.structure_density*constants.greenhouse.vol_cap;

constants.greenhouse.structure_total_mass = constants.greenhouse.mass_cap + constants.greenhouse.mass_inner + constants.greenhouse.mass_outer;

constants.greenhouse.VolWorking = 3.1415*constants.greenhouse.radius_inner^2 * constants.greenhouse.length;

constants.greenhouse.AreaWorking = 2*3.1415*constants.greenhouse.radius_avg^2 + 2*3.1415*constants.greenhouse.radius_avg*constants.greenhouse.length;

constants.greenhouse.atm_mass = constants.greenhouse.atm_density * constants.greenhouse.VolWorking;

constants.greenhouse.HeatMin = determineinternalHeat(constants.greenhouse.VolWorking, constants.greenhouse.min_temp, constants);
constants.greenhouse.HeatMax = determineinternalHeat(constants.greenhouse.VolWorking, constants.greenhouse.max_temp, constants);

constants.greenhouse.lightWatts = 32; % Watts. How much power the lights are using, current estimate is based on 2 32 Watt bulbs.

constants.greenhouse.HeatGen = constants.greenhouse.lightWatts;

constants.greenhouse.thickness_max = 0.02; 


return;

end