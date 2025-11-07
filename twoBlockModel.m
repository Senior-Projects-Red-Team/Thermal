function [Q_dots, heating, powerUse] = twoBlockModel(Qs,Masses, Cps, Rs, Q_gen,Q_Lights, heatersMax, t)
Qlights_Actual = Q_Lights;
% TODO: Add light cycle
Q_dots = zeros(length(Qs),1);

%% Unpacking Heaters
Rs_Heaters = Rs(1:4);
Cps_Heaters = Cps(1:4);
Qs_Heaters = Qs(1:4);
Masses_Heaters = Masses(1:4);
Ts_Heaters = Qs_Heaters./(Masses_heaters.*Cps_Heaters);

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
Qdots_Structure = zeros(4:1);

%% Unpacking Heater Destinations
% Water = 11 (Note here R is from the water to the Life Support Module) -> 1
% Electronics Block = 12 -> 2
% Life Support Block = 13 -> 3

Rs_Destinations = Rs(11:13);
Cps_Destinations = Cps(11:13);
Qs_Destinations = Qs(11:13);
Masses_Destinations = Masses(11:13);
Ts_Destinations = Qs_Destinations./(Masses_Destinations.*Cps_Destinations);\
Qdots_Destinations = zeros(3:1);

%% Unpacking Moon
Rs_Moon = Rs(14:end);
Cps_Moon = Cps(14:end);
Qs_Moon = Qs(14:end);
Masses_Moon = Masses(14:end);
Ts_Moon = Qs_Moon./(Masses_Moon.*Cps_Moon);
Qdots_Moon = zeros(length(Rs_Moon),1);



%% Lighting Effect

%% Resistance Network

% Internals (Destinations and Plate)
 % Heaters
Qdots_Destinations(1:2) = Qdots_Destinations(1:2) + (Ts_Heaters(1:2) - Ts_Destinations(1:2))./(Rs_Heaters(1:2));
Qdots_Destinations(3) = Qdots_Destinations(3) + (Ts_Heaters(3) - Ts_Destinations(3))./Rs_Heaters(3);
Qdots_Destinations(3) = Qdots_Destinations(3) + (Ts_Heaters(4) - Ts_Destinations(3))./Rs_Heaters(4);
 % Water to Life Support (and Visa Versa)
Qdots_Destinations(1) = Qdots_Destinations(1) - (Ts_Destinations(1) - Ts_Destinations(3))./Rs_Destinations(1);
Qdots_Destinations(3) = Qdots_Destinations(3) + (Ts_Destinations(1) - Ts_Destinations(3))./Rs_Destinations(1);

 % Eletcronics Block
Qdots_Destinations(2) = Qdots_Destinations(2) + Qgen; % Effect of Electronics
Qdots_Destinations(2) = Qdots_Destinations(2) - (Ts_Destinations(2) - Ts_Structure(2))./Rs_convWall; % Convection To Walls
Qdots_Structure(2) = Qdots_Structure(2) + (Ts_Destinations(2) - Ts_Structure(2))./Rs_convWall;
Qdots_Destinations(2) = Qdots_Destinations(2) - (Ts_Destinations(2) - Ts_Structure(1))./Rs_convPlate; % Convection To Plate
Qdots_Structure(1) = Qdots_Structure(2) + (Ts_Destinations(2) - Ts_Structure(1))./Rs_convWall;

 % Life Support Module
Qdots_Destinations(3) = Qdots_Destinations(3) + Qlights_Actual; % Effect of Lights
Qdots_Destinations(3) = Qdots_Destinations(3) - (Ts_Destinations(3) - Ts_Structure(2))./Rs_convWall; % Convection To Walls
Qdots_Structure(2) = Qdots_Structure(2) + (Ts_Destinations(3) - Ts_Structure(2))./Rs_convWall;
Qdots_Destinations(3) = Qdots_Destinations(3) - (Ts_Destinations(3) - Ts_Structure(1))./Rs_convPlate; % Convection To Plate
Qdots_Structure(1) = Qdots_Structure(2) + (Ts_Destinations(3) - Ts_Structure(1))./Rs_convWall;

 

% Externals (Through the Wall and into the moon)
 % Wall
Qdots_Structure(2) = Qdots_Structure(2) - (Ts_Structure(2) - Ts_Structure(3))/Rs_Structure(2); % Inner Wall to Bridges
Qdots_Structure(3) = Qdots_Structure(3) + (Ts_Structure(2) - Ts_Structure(3))/Rs_Structure(2) - (Ts_Structure(3) - Ts_Structure(4))/Rs_Structure(3); % Bridges
Qdots_Structure(4) = Qdots_Structure(4) + (Ts_Structure(3) - Ts_Structure(4))/Rs_Structure(3) - (Ts_Structure(4) - Ts_Moon(1))/Rs_Structure(4); % Bridges

 % Moon
Qdots_Moon(1) = Qdots_Moon(1) + (Ts_Structure(4) - Ts_Moon(1))/Rs_Structure(4) - (Ts_Moon(1) - Ts_Moon(2))./Rs_Moon(1);
for i = 2:length(Rs_Moon)-1
    Qdots_Moon(i) = Qdots_Moon(i) + (Ts_Moon(i-1) - Ts_Moon(i))/Rs_Moon(i-1) - (Ts_Moon(i) - Ts_Moon(i+1))./Rs_Moon(i);
end
Qdots_Moon(end) = Q_dots_Moon(end) + (Ts_Moon(end-1) - Ts_Moon(end))/Rs_Moon(end-1);


% PD Control
heating = zeros(4,1);
for i = 1:3
    if (Ts_Destinations(i) < 25+273.15)
        deltaT = Ts_Destinations(i) - (25+273.15);
        reldT = deltaT / 3;
        relQdots = Qdots_Destinations(i) / (2*heatersMax);
        heating(i) = heating(i) + heatersMax * (relQdots+reldT);
        if(i == 3) % As we have 2 heaters for this destination we treat heatersMax as 2x its value and apply to both heaters 3 and 4
            relQdots = Qdots_Destinations(i) / (4*heatersMax);
            heating(3:4) = heating(3:4) + heatersMax * (relQdots+reldT);
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
Q_dots(14:end) = Qdots_Moon;


end