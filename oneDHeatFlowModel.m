function [Q_dots] = oneDHeatFlowModel(Qs,Masses, Cps, Rs, Q_gen, t)
%This function finds the rate of change of heat in each slice, for use in
%ODE45

Temps = Qs./(Masses.*Cps);

Entries = length(Qs);

Q_dots = zeros(length(Qs),1);

Q_dots(1) = Q_gen - (Temps(1)-Temps(2))/Rs(1);

for i= 2:(Entries-1)
Q_dots(i) = (Temps(i-1)-Temps(i))/Rs(i-1) - (Temps(i)-Temps(i+1))/Rs(i);
end

Q_dots(Entries) = (Temps(Entries-1)-Temps(Entries))/Rs(Entries-1);

Q_dots = Q_dots./1000;


end