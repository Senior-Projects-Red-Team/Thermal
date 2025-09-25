function [Q_dots] = oneDWithHeatPump(Qs,Masses, Cps, Rs, Q_gen, pump_Max, heatersBool, t)
%This function finds the rate of change of heat in each slice, for use in
%ODE45

Temps = Qs./(Masses.*Cps);

heating_max = 0.9*Q_gen;

if( mod((t/3600),24) > 16 )
    Q_gen = 0.1* Q_gen;
end


Entries = length(Qs);

Q_dots = zeros(length(Qs),1);

Q_dots(1) = Q_gen - (Temps(1)-Temps(2))/Rs(1);



if(Temps(1) < 24 + 273.15 && Q_dots(1) < 0 && heatersBool)
    heating = -1 * Q_dots(1);
    if(heating > heating_max)
        heating = heating_max;
    end
    Q_dots(1) = Q_dots(1) + heating;
end



if(Temps(1) > (26 + 273.15) && Q_dots(1) > 0 && pump_Max ~= 0)
    to_Pump = Q_dots(1);

    if(to_Pump >= pump_Max)
        to_Pump = pump_Max;
    end

    if(to_Pump < 0)
        to_Pump = 0;
    end

    Q_dots(1) = Q_dots(1) - to_Pump;
    Q_dots(4) = Q_dots(4) + to_Pump;
end


for i= 2:(Entries-1)
Q_dots(i) = Q_dots(i) + (Temps(i-1)-Temps(i))/Rs(i-1) - (Temps(i)-Temps(i+1))/Rs(i);
end

Q_dots(Entries) = (Temps(Entries-1)-Temps(Entries))/Rs(Entries-1);




Q_dots = Q_dots./1000; % Converts from W to kW



end
