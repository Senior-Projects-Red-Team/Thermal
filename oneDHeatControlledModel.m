function [ts_actual, Qs_actual, Q_dots, heating, powerUse] = oneDHeatControlledModel(q_gen_base, q_gen_lights, powerMax, heatPumpMax, heatersOn, network)

masses = network.masses;
cps = network.cps;
Rs = network.Rs;
Heats = network.Heats;

Q_0 = Heats;
ts_actual = zeros(4e6,1);
Qs_actual = zeros(4e6,10);
t_i = 1;
for i = 1:15*24
[ts,Qs] = ode45(@(t,Qs) oneDWithHeatPump(Qs, masses, cps, Rs, q_gen_base, q_gen_lights, heatPumpMax, powerMax, heatersOn, t),[(i-1)*3600,i*3600],Q_0);
ts_actual(t_i: length(ts) + (t_i - 1)) = ts;
Qs_actual(t_i: length(ts) + (t_i - 1),:) = Qs(:,1:10);
t_i = t_i + length(ts);
Q_0 = Qs(end,:);
end
toc()

ts_actual = ts_actual(1:50:t_i-1);
Qs_actual = Qs_actual(1:50:t_i-1,:);

%Extract extra variables by re-evaluating the ODE function
Q_dots = 0.*Qs_actual;
heating = zeros(length(ts_actual),1);
powerUse = zeros(length(ts_actual),1);
for i = 1:length(ts_actual)
    [Q_dots(i,:), heating(i), powerUse(i)] = oneDWithHeatPump(transpose(Qs_actual(i,:)), masses(1:10), cps(1:10), Rs(1:10), q_gen_base, q_gen_lights, heatPumpMax, powerMax, heatersOn, ts_actual(i));
    % floor(100*i/length(ts_actual))
    % toc()
end
toc()

end

