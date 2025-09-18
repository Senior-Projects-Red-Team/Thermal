function constants = loadConstants()

constants.regolith.conductivity = 2.95e-2; % 2.95e-2 Watts/m*K (Worst case scenario, revise as needed) % Maybe 7.4 e-4 W/mK @SP Surface 3.4e-3 W/mK @ 1m ??? % Range 0.5e-3 -> 0.8e-3

constants.regolith.Mintemp = 40; % K
constants.regolith.Maxtemp = 125; % K
constants.regolith.density = 1800; % kg/m^3 @ 1m
constants.regolith.cp = 840; % J/kgK
constants.regolith.diffusivity = constants.regolith.conductivity / (constants.regolith.density * constants.regolith.cp); % m^2/s; Seems way too low.
constants.regolith.diffusivity_alt = 2e-7; % m^2/s
constants.regolith.conductivity_alt = constants.regolith.diffusivity_alt*constants.regolith.density*constants.regolith.cp;% W/m*K

constants.power.battery_mass = 0;
constants.power.battery_cp = 0.8; % kj/kgK 

constants.greenhouse.max_temp = 28 + 273.15; % K
constants.greenhouse.min_temp = 22 + 273.15; % K
constants.greenhouse.atm_cp = 1.005; % kJ/kgK % Assuming similar to normal air for now
constants.greenhouse.structure_cp = 0.921 ; % 0.921 for Aluminum, 0.75 for Carbon Fiber, 1.5 for Kevlar, HDPE 1.3-2.2, all units as kJ/kgK
constants.greenhouse.water_cp = 4.18; % kJ/kgK
constants.greenhouse.atm_density = 1.204; % kg/m^3 Assumes normal air
constants.greenhouse.water_density = 1000; %kg/m^3
constants.greenhouse.water_mass = 1; % I'm temporarily assuming 1 kg of water, more specific values needed.
constants.greenhouse.structure_mass = 6; % Consult with SE and Structures Teams
constants.greenhouse.AreaMin = 2*0.2*0.2 + 4*0.2*0.1; % 20cm tall, 20cm long, 10cm wide %% Final answer in m^2
constants.greenhouse.AreaMax = 2*0.5*0.5 + 4*0.5*0.25; % 50 cm, 50cm, 25cm %% Final answer in m^2
constants.greenhouse.VolMin = 0.2*0.2*0.1; % m^3
constants.greenhouse.VolMax = 0.5*0.5*0.25; % m^3

constants.greenhouse.AreaWorking = 6*0.35*0.35;
constants.greenhouse.VolWorking = 0.35^3;

constants.greenhouse.HeatMin = determineinternalHeat(constants.greenhouse.VolWorking, constants.greenhouse.min_temp, constants);
constants.greenhouse.HeatMax = determineinternalHeat(constants.greenhouse.VolWorking, constants.greenhouse.max_temp, constants);

constants.greenhouse.lightWatts = 32; % Watts. How much power the lights are using, current estimate is based on 2 32 Watt bulbs.

constants.greenhouse.HeatGen = constants.greenhouse.lightWatts;

constants.greenhouse.thickness_max = 0.02; 


return;

end