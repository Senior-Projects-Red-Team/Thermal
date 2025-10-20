function runModelTemps(Masses, Cps, Rs, Heats, Q_gen, heatPump, heatersBool, wallBool, surfaceBool, slicesBool, figNum)

figure(figNum)
hold on
for i = 1:(15*24)

    if(i == 1)
        [ts,Qs] = ode45(@(t,Qs) oneDWithHeatPump(Qs, Masses, Cps, Rs, Q_gen, heatPump, heatersBool, t),[0,i*3600],Heats);
    else
        [ts,Qs] = ode45(@(t,Qs) oneDWithHeatPump(Qs, Masses, Cps, Rs, Q_gen, heatPump, heatersBool, t),[(i-1)*3600,i*3600],Qs(end,:));
    end
    plot(ts./(24*3600),Qs(:,1)./(Masses(1)*Cps(1))-273.15,Color="b")

    if(wallBool)
        for j = 2:4
            plot(ts./(24*3600),Qs(:,j)./(Masses(j)*Cps(j))-273.15,Color="r",LineStyle="-.")
        end
    end

    if(surfaceBool)
        plot(ts./(24*3600),Qs(:,5)./(Masses(5)*Cps(5))-273.15,Color="r",LineStyle="-.")
    end
    if(slicesBool)
        for j = 5:length(Rs)
            if(mod(j,5)==0)
                plot(ts./(24*3600),Qs(:,j)./(Masses(j)*Cps(j))-273.15,LineStyle=":",Color="g")
            end
        end
    end

    if(i == 1)
        yline(28)
        yline(22)
        title("Temp of Internal Atmosphere with " + Q_gen + " (W) Internal Heat Gen and" + heatPump + "(W) Heat Pump")
        ylabel("C")
        xlabel("Time (Days)")
    end

end
legend(["Internal", "Lunar Surface", "Viable Range"],Location="southeast")
hold off

end
