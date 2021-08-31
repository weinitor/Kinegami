% Creating a Crease Schematic - Origami elbow fitting
% Last edited 6/8/2021 by Lucien Peach

clear
close all
clc

% Specify inputs
r = 0.02; %[m]
n = 4;
<<<<<<< HEAD
phi = pi/2; %[rad]
=======
phi = 3.14/3; %[rad]
>>>>>>> 3a8c5293a489e141dc86328bef9f1d0e1f607749
theta = 3*pi/2; %[rad] - splits for greater than pi/2

% Automatically @pi/2 split unless specified by user to not split
split = 'off';
mirror = 'on';

% Outputs array of lenghts and value of ls in [m]
<<<<<<< HEAD
[lengths, ls] = Origami_Elbow_Parameters(r, n, phi, theta, split);
=======
[lengths, ls] = Origami_Elbow_Parameters(r, n, phi, theta);
TuckAngles(r, n, phi, theta)
>>>>>>> 3a8c5293a489e141dc86328bef9f1d0e1f607749

% Create a figure that demonstrates the crease schematic

% Specify values for h1 and h2, the heights of the two tube sections
h1 = 0.1; %[m]
h2 = 0.1; %[m]

% Outputs graphing for elbow fitting
[dataFoldA, m, lmax] = Origami_Elbow_CreasePattern(lengths, ls, n, h1, h2, r, phi, theta, mirror, split);
axis off

% Convert to DXF
% filename = (['FoldA_r' num2str(r) '_n' num2str(n) '_phi' num2str(phi) ...
%     '_theta' num2str(theta) '.dxf']);
% GenerateDXF(filename, dataFoldA)
