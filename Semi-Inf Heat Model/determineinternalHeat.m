function heat = determineinternalHeat(volume, temperature, constants)
% Determines the minimum expected internal heat energy properties of the
% Greenhouse.
%
heat = volume * temperature * constants.greenhouse.atm_cp; % accounts for air
heat = heat + constants.greenhouse.water_mass*temperature*constants.greenhouse.water_cp; % Accounts for water
heat = heat + constants.greenhouse.structure_mass*temperature*constants.greenhouse.structure_cp; % Accounts for the structure



end