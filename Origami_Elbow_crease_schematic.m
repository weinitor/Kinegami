% Creating a Crease Schematic - Fold Type A
% Last edited 6/8/2021 by Lucien Peach

clear
clc

% Specify inputs
r = 0.02; %[m]
n = 4;
phi = 3.14; %[rad]
theta = pi; %[rad] - splits for greater than pi/2
mirror = 'on';

% Outputs array of lenghts and value of ls in [m]
[lengths, ls] = Origami_Elbow_creasedesign(r, n, phi, theta);

% Create a figure that demonstrates the crease schematic

% Specify values for h1 and h2, the heights of the two tube sections
h1 = 0.1; %[m]
h2 = 0.1; %[m]

% Outputs graphing for elbow fitting
[dataFoldA, m, lmax] = Origami_Elbow_papercut(lengths, ls, n, h1, h2, r, phi, theta, mirror);
axis off

% Convert to DXF
% filename = (['FoldA_r' num2str(r) '_n' num2str(n) '_phi' num2str(phi) ...
%     '_theta' num2str(theta) '.dxf']);
% GenerateDXF(filename, dataFoldA)
