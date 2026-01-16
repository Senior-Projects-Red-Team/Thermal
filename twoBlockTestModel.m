function [Q_dots, heating, powerUse] = twoBlockTestModel(Qs, Masses, Cps, Rs, Q_gen, Q_Lights, heatersMax, T_ext , t)
Qlights_Actual = Q_Lights; % Catches Daytime
% TODO: Add light cycle

sunset = 0.1; % Length of a sunset/rise in minutes

t_eff = mod(t,(24*3600)); % gets time in 24 hour increments
t_eff_hours = t_eff/3600; % Used for debugging and integral control.

if(t_eff >= 16*3600) % Catches Sunset
    light_percentage = 1-(t_eff - 16*3600)/(60*sunset);
    Qlights_Actual = Q_Lights * light_percentage;
end

if(t_eff >= 16*3600 + 60*sunset) % Post-Sunset
    Qlights_Actual = 0;
end

if(t_eff >= 24*3600 - 60 *sunset) % Catches Sunrise
    light_percentage = (t_eff - (24*3600 - 60 *sunset))/(60*sunset);
    Qlights_Actual = Q_Lights * light_percentage;
end

% if( mod((t/3600),24) > 16 ) % Turns off the lights for 8 hours a day.
%     Qlights_Actual = 0;
% end

Q_dots = zeros(length(Qs),1);

%% Unpacking Heaters
Rs_Heaters = Rs(1:4);
Cps_Heaters = Cps(1:4);
Qs_Heaters = Qs(1:4);
Masses_Heaters = Masses(1:4);
Ts_Heaters = Qs_Heaters./(Masses_Heaters.*Cps_Heaters);

%% Unpacking Structure
% Plate = 5; -> 1
% innerWall = 6; -> 2
% Bridges = 7; -> 3
% outerWall = 8; -> 4
% Conv to Plate = 9 % No mass props
% Conv to Walls = 10 % No mass props
Rs_Structure = Rs(5:8);
Rs_convPlate = Rs(9);
Rs_convWall = Rs(10);
Cps_Structure = Cps(5:8);
Qs_Structure = Qs(5:8);
Masses_Structure = Masses(5:8);
Ts_Structure = Qs_Structure./(Masses_Structure.*Cps_Structure);
Qdots_Structure = zeros(4,1);

%% Unpacking Heater Destinations
% Water = 11 (Note here R is from the water to the Life Support Module) -> 1
% Electronics Block = 12 -> 2
% Life Support Block = 13 -> 3

Rs_Destinations = Rs(11:13);
Cps_Destinations = Cps(11:13);
Qs_Destinations = Qs(11:13);
Masses_Destinations = Masses(11:13);
Ts_Destinations = Qs_Destinations./(Masses_Destinations.*Cps_Destinations);
Qdots_Destinations = zeros(3,1);



%% Lighting Effect

%% Resistance Network

% Internals (Destinations and Plate)
 % Heaters
Qdots_Destinations(1:2) = Qdots_Destinations(1:2) + (Ts_Heaters(1:2) - Ts_Destinations(1:2))./(Rs_Heaters(1:2)); % This step models heat flow from heaters to the water and electronics modules
Qdots_Destinations(3) = Qdots_Destinations(3) + (Ts_Heaters(3) - Ts_Destinations(3))./Rs_Heaters(3); % These two steps model heat flow from heaters to the life support module
Qdots_Destinations(3) = Qdots_Destinations(3) + (Ts_Heaters(4) - Ts_Destinations(3))./Rs_Heaters(4);
 % Water to Life Support (and Visa Versa)
Qdots_Destinations(1) = Qdots_Destinations(1) - (Ts_Destinations(1) - Ts_Destinations(3))./Rs_Destinations(1);
Qdots_Destinations(3) = Qdots_Destinations(3) + (Ts_Destinations(1) - Ts_Destinations(3))./Rs_Destinations(1);

 % Eletcronics Block
Qdots_Destinations(2) = Qdots_Destinations(2) + Q_gen; % Effect of Electronics
Qdots_Destinations(2) = Qdots_Destinations(2) - (Ts_Destinations(2) - Ts_Structure(2))./Rs_convWall; % Convection To Walls
Qdots_Structure(2) = Qdots_Structure(2) + (Ts_Destinations(2) - Ts_Structure(2))./Rs_convWall;
Qdots_Destinations(2) = Qdots_Destinations(2) - (Ts_Destinations(2) - Ts_Structure(1))./Rs_convPlate; % Convection To Plate
Qdots_Structure(1) = Qdots_Structure(1) + (Ts_Destinations(2) - Ts_Structure(1))./Rs_convPlate;

 % Life Support Module
Qdots_Destinations(3) = Qdots_Destinations(3) + Qlights_Actual; % Effect of Lights
Qdots_Destinations(3) = Qdots_Destinations(3) - (Ts_Destinations(3) - Ts_Structure(2))./Rs_convWall; % Convection To Walls
Qdots_Structure(2) = Qdots_Structure(2) + (Ts_Destinations(3) - Ts_Structure(2))./Rs_convWall;
Qdots_Destinations(3) = Qdots_Destinations(3) - (Ts_Destinations(3) - Ts_Structure(1))./Rs_convPlate; % Convection To Plate
Qdots_Structure(1) = Qdots_Structure(1) + (Ts_Destinations(3) - Ts_Structure(1))./Rs_convPlate;

 

% Externals (Through the Wall and into the moon)
 % Wall
Qdots_Structure(2) = Qdots_Structure(2) - (Ts_Structure(2) - Ts_Structure(3))/Rs_Structure(2); % Inner Wall to Bridges
Qdots_Structure(3) = Qdots_Structure(3) + (Ts_Structure(2) - Ts_Structure(3))/Rs_Structure(2) - (Ts_Structure(3) - Ts_Structure(4))/Rs_Structure(3); % Bridges
Qdots_Structure(4) = Qdots_Structure(4) + (Ts_Structure(3) - Ts_Structure(4))/Rs_Structure(3) - (Ts_Structure(4) - (T_ext))/Rs_Structure(4); % Outer Wall


% PD Control

T_maintain = 25 + 273.15;

heating = zeros(4,1);
for i = 1:3
    if (Ts_Destinations(i) < T_maintain)
        deltaT = (T_maintain)-Ts_Destinations(i);
        reldT = deltaT / 3;
        rel_Int = 0;
        %rel_Int = (((t/(24*3600))+2) * reldT); % Implements integral control. Comment out when not desired.
        relQdots = Qdots_Destinations(i) / (2*heatersMax);
        heating(i) = heating(i) + heatersMax * (relQdots+reldT + rel_Int);
        if(i == 3)
            relQdots = Qdots_Destinations(3) / (4*heatersMax);
            heating(3) = heatersMax * (relQdots+reldT + rel_Int);
            heating(4) = heatersMax * (relQdots+reldT + rel_Int);
        end
    end
end
for i = 1:4
    if(heating(i) < 0)
        heating(i) = 0;
    end
    if(heating(i) > heatersMax)
        heating(i) = heatersMax;
    end
end
% Heaters
Qdots_Heaters = heating(1:4);
Qdots_Heaters(1:2) = Qdots_Heaters(1:2) - (Ts_Heaters(1:2) - Ts_Destinations(1:2))./(Rs_Heaters(1:2));
Qdots_Heaters(3) = Qdots_Heaters(3) - (Ts_Heaters(3) - Ts_Destinations(3))./Rs_Heaters(3);
Qdots_Heaters(4) = Qdots_Heaters(4) - (Ts_Heaters(4) - Ts_Destinations(3))./Rs_Heaters(4);

powerUse = sum(heating) + Q_gen + Qlights_Actual;

% Repacking
Q_dots(1:4) = Qdots_Heaters;
Q_dots(5:8) = Qdots_Structure;
Q_dots(9:10) = 0; % Handles blank space for Convections
Q_dots(11:13) = Qdots_Destinations;

Q_dots = Q_dots/1000; % Converts from W to kW

end