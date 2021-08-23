% Creating a Crease Schematic - Origami Revolute Joint
% Last edited 6/17/2021 by Lucien Peach

clear
close all
clc

% Specify inputs (3D Modeling will also consider w, but do not worry about
% this value for the time being for 2D)
r = 0.02; %[m]
n = 6; % must be even, 4 or greater
theta_m = (3/2)*pi; %[rad]

% Outputs array of lengths and value of ls in [m]
[lengths, ls] = Origami_RevoluteJoint_Parameters(r, n, theta_m);

% Create a figure that demonstrates the crease schematic

% Specify values for h1 and h2, the heights of the two tube sections
h1 = 0.03; %[m]
h2 = 0.03; %[m]
nz = 4;

% Outputs graphing for elbow fitting
[dataFoldD, m, lmax] = Origami_RevoluteJoint_CreasePattern(lengths, ls, n, h1, h2, r, theta_m, nz);
axis off

% Convert to DXF
% filename = (['Revolute_r' num2str(r) '_n' num2str(n) '_theta' num2str(theta_m) ...
%     '.dxf']);
% GenerateDXF(filename, dataFoldD)
