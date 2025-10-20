function runModelHeats(Masses, Cps, Rs, Heats, Q_gen, figNum, constants)

fig = figNum+10;

figure(fig)
hold on
for i = 1:(15*24)

    if(i == 1)
        [ts,Qs] = ode45(@(t,Qs) oneDHeatFlowModel(Qs, Masses, Cps, Rs, Q_gen, t),[0,i*3600],Heats);
    else
        [ts,Qs] = ode45(@(t,Qs) oneDHeatFlowModel(Qs, Masses, Cps, Rs, Q_gen, t),[(i-1)*3600,i*3600],Qs(end,:));
    end

    greenhouse_heat = Qs(:,1) + Qs(:,2) + Qs(:,3) + Qs(:,4); 

    plot(ts./(24*3600),greenhouse_heat,Color="b")
end

yline(constants.greenhouse.HeatMin)
yline(constants.greenhouse.HeatMax)
title("Total Greenhouse Heat over time with " + Q_gen + " (W) Internal Heat Gen")
ylabel("Heat (kJ)")
xlabel("Time (Days)")

hold off

end