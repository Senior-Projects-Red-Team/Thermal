clc;
clear;
close all;

% Units in meters
r_0 = 0.17;
pi = 3.1415;
roll_thickness = 0.005;

sum = 0;

while r_0 <= 0.17 + 0.032
    sum = sum + 2*pi*r_0;
    r_0 = r_0 + roll_thickness;
end


