% ORIGAMI_TWIST_STANDALONE - Generates the crease pattern of the origami
% twist fitting. 

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last Edited 12/15/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


close all
clear 
clc

% Add all the folders and subfolders to the search path
addpath(genpath(fileparts(mfilename('fullpath'))));

% Specify inputs
r = 0.02; %[m]
n = 6;
h = 0.03; %[m]
alpha = pi/4; %[rad]

% Outputs values for x and l in [m]
[x, l, ls] = Origami_Twist_Parameters(r, n, h, alpha);

% Create a figure that demonstrates the crease schematic

% Specify values for h1 and h2, the heights of the two tube sections
h1 = 0.03; %[m]
h2 = 0.03; %[m]

% Outputs midsection for graphing
[dataFoldB, m, lmax] = Origami_Twist_CreasePattern(x, l, ls, n, h1, h2, alpha);

% Create Duplication for Overlap Slide
[dataFoldNew] = StandaloneDuplication(dataFoldB, ls, n, lmax, 'twist', h1);

% Plotting
figure()
hold on
for i = 1:size(dataFoldNew, 2)
    
    plot(dataFoldNew(i).x, dataFoldNew(i).y, 'color', dataFoldNew(i).color)    
    
end

% Label the plot for clarity
title({
    ('Twist Joint for Provided Parameters:')
    ['[r = ' num2str(r) ', n = ' num2str(n) ', h = ' num2str(h) ', alpha = ' num2str(alpha) ']']
    })

daspect([1 1 1])
axis off
set(gcf, 'color', 'w')

% Convert to DXF
filename = (['Twist_r' num2str(r) '_n' num2str(n) '_h' num2str(h) ...
    '_alpha' num2str(alpha) '.dxf']);
GenerateDXF(filename, dataFoldNew)
