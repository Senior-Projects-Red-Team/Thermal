
% Created 9/5/2025

clear;
close all;
clc;

%% Load Constants
tic();
constants = loadConstants();





%% Temp of external surface

Q_dot = constants.greenhouse.HeatGen; % Here we're assuming that we'll maintain the internal temperature.

%semiInfModel(32, constants); % Low end
semiInfModel(50, constants); % Theoretical maximum given by Power Team
%semiInfModel(64, constants);
%semiInfModel(80, constants);
%semiInfModel(96, constants);
%semiInfModel(128, constants); % Approximates Earth Sunlight.


xs = linspace(0,3,601);

T_ratio = erfc(xs./sqrt(4*constants.regolith.diffusivity_alt*(15*24*60*60)));

plot(xs, T_ratio)


toc()

