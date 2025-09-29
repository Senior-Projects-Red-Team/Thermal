function [Q_dots] = oneDWithHeatPump(Qs,Masses, Cps, Rs, Q_gen,Q_Lights, pump_Max, powerMax, heatersBool, t)
%This function finds the rate of change of heat in each slice, for use in
%ODE45

Temps = Qs./(Masses.*Cps);

Q_gen_total = Q_gen + Q_Lights;

if( mod((t/3600),24) > 16 )
    Q_gen_total = Q_gen;
end

Entries = length(Qs);

Q_dots = zeros(length(Qs),1);

Q_dots(1) = Q_gen_total - (Temps(1)-Temps(2))/Rs(1);


if(Temps(1) < 25 + 273.15 && Q_dots(1) < 0 && heatersBool)
    heatMax = powerMax - Q_gen_total;
    damping_factor = -Q_dots(1)/(2*heatMax);
    temp_offset = abs((25 + 273.15) - Temps(1)); % Captures the amount below 25 C we are
    offset_factor = temp_offset / 3; % Normalizes the offset from 0 -> 1
    heating_factor = offset_factor + damping_factor;
    if(heating_factor > 1)
        heating_factor = 1;
    end
    heating = heating_factor * heatMax;

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
