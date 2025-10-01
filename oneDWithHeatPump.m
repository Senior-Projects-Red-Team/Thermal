function [Q_dots, heating, powerUse] = oneDWithHeatPump(Qs,Masses, Cps, Rs, Q_gen,Q_Lights, pump_Max, powerMax, heatersBool, t)
%This function finds the rate of change of heat in each slice, for use in
%ODE45

heating = 0;

Temps = Qs./(Masses.*Cps);

Q_gen_total = Q_gen + Q_Lights;

if( mod((t/3600),24) > 16 )
    Q_gen_total = Q_gen;
end

Entries = length(Qs);

Q_dots = zeros(length(Qs),1);

Q_dots(1) = Q_gen_total - (Temps(1)-Temps(2))/Rs(1);


if(Temps(1) < 25 + 273.15 && heatersBool)
    heatMax = powerMax - Q_gen_total;
    tempDiff = (25+273.15) - Temps(1);
    tempDiffRatio = tempDiff / 2.5; % Normalized Offset factor (k)
    Q_dot_ratio = -1*Q_dots(1) / (3*heatMax); % Normalized Damping factor (d) % Creates noise if too high
    
    heating = heatMax*(tempDiffRatio + Q_dot_ratio);

    if(heating > heatMax)
        heating = heatMax;
    end
    if(heating < 0)
        heating = 0;
    end

end

Q_dots(1) = Q_dots(1) + heating;

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

powerUse = heating + Q_gen_total;

if(mod(floor(t), 2*3600) == 0)
    breakpoint = 1;
end

Q_dots = Q_dots./1000; % Converts from W to kW



end
