% ORIGAMI_FINGERTIP_STANDALONE - Generates the crease pattern of the 
% origami fingertip "joint". 

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last Edited 8/3/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


clear
close all
clc

% Add all the folders and subfolders to the search path
addpath(genpath(fileparts(mfilename('fullpath'))));

% Specify inputs (3D Modeling will also consider w, but do not worry about
% this value for the time being for 2D)
r = 0.02; %[m]
n = 4; % must be even, 4 or greater
theta_m = (3/2)*pi; %[rad]

% Outputs array of lengths and value of ls in [m]
[lengths, ls] = Origami_RevoluteJoint_Parameters(r, n, theta_m);

% Create a figure that demonstrates the crease schematic

% Specify value for h1. 
h1 = 0; %[m]

% Outputs graphing for elbow fitting
[dataFoldF, m, lmax] = Origami_Fingertip_CreasePattern(lengths, ls, n, h1);

% Create Duplication for Overlap Slide
[dataFoldNew] = StandaloneDuplication(dataFoldF, ls, n, lmax, 'fingertip', h1);

% Plotting
figure()
hold on
for i = 1:size(dataFoldNew, 2)
    
    plot(dataFoldNew(i).x, dataFoldNew(i).y, 'color', dataFoldNew(i).color)    
    
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
filename = (['Revolute_r' num2str(r) '_n' num2str(n) '_theta' num2str(theta_m) ...
    '.dxf']);
GenerateDXF(filename, dataFoldNew)
