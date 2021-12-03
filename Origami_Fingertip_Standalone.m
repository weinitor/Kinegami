% Creating a Crease Schematic - Fingertip "Joint"
% Last edited 8/3/2021 by Lucien Peach

clear
close all
clc

% Specify inputs (3D Modeling will also consider w, but do not worry about
% this value for the time being for 2D)
r = 0.02; %[m]
n = 8; % must be even, 4 or greater
theta_m = (3/2)*pi; %[rad]

% Outputs array of lengths and value of ls in [m]
[lengths, ls] = Origami_RevoluteJoint_Parameters(r, n, theta_m);

% Create a figure that demonstrates the crease schematic

% Specify values for h1 and h2, the heights of the two tube sections
h1 = 0.03; %[m]
nz = 1;

% Outputs graphing for elbow fitting
[dataFoldF, m, lmax] = Origami_Fingertip_CreasePattern(lengths, ls, n, h1, r, theta_m);

% Plotting
figure()
hold on
for i = 1:size(dataFoldF, 2)
    
    plot(dataFoldF(i).x, dataFoldF(i).y, 'color', dataFoldF(i).color)    
    
end

% Label the plot for clarity
title({
    ('Origami Schematic 2.A for Provided Parameters:')
    ['[r = ' num2str(r) ', n = ' num2str(n) ', theta = ', num2str(theta_m) ']']
    })

daspect([1 1 1])
axis off
set(gcf, 'color', 'w')

% Convert to DXF
% filename = (['Revolute_r' num2str(r) '_n' num2str(n) '_theta' num2str(theta_m) ...
%     '.dxf']);
% GenerateDXF(filename, dataFoldD)
