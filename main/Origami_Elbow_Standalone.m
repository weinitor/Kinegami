% ORIGAMI_ELBOW_STANDALONE - Generates the crease pattern of the origami
% elbow fitting. 

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last Edited 6/8/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


clear
close all
clc

% Add all the folders and subfolders to the search path
addpath(genpath(fileparts(mfilename('fullpath'))));

% Specify inputs
r = 0.02; %[m]
n = 4;

phi = pi/3; %[rad]
theta = 2*pi/3; %[rad] - splits for greater than pi/2

% phi = pi/2; %[rad]
% theta = -3*pi/4; %[rad] - splits for greater than pi/2

% Automatically @pi/2 split unless specified by user to not split
split = 'on';
tuck = 'on';

% Outputs array of lengths and value of ls in [m]
[lengths, ls] = Origami_Elbow_Parameters(r, n, phi, theta, split);
[tuckangle] = TuckAngles(r, n, phi, theta, split);

% Create a figure that demonstrates the crease schematic

% Specify values for h1 and h2, the heights of the two tube sections
h1 = 0.02; %[m]
h2 = 0.02; %[m]

% Outputs graphing for elbow fitting
[dataFoldA, m, lmax] = Origami_Elbow_CreasePattern(lengths, ls, n, h1, h2, ...
    theta, tuck, split, tuckangle);

% Create Duplication for Overlap Slide
[dataFoldNew] = StandaloneDuplication(dataFoldA, ls, n, lmax, 'elbow', h1);

% Plotting
figure()
hold on
for i = 1:size(dataFoldNew, 2)
    
    plot(dataFoldNew(i).x, dataFoldNew(i).y, 'color', dataFoldNew(i).color)    
    
end

% Label the plot for clarity
title({
    ('Elbow Joint for Provided Parameters:')
    ['[r = ' num2str(r) ', n = ' num2str(n) ', phi = ' num2str(phi) ', theta = ' num2str(theta) ']']
    })

daspect([1 1 1])
axis off
set(gcf, 'color', 'w')

% Convert to DXF
filename = (['Elbow_r' num2str(r) '_n' num2str(n) '_phi' num2str(phi) ...
    '_theta' num2str(theta) '.dxf']);
GenerateDXF(filename, dataFoldNew)

