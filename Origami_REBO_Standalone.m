% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.txt for detail.
% Authors: 
% Lucien Peach <peach@seas.upenn.edu>

% Creating a Crease Schematic - Origami REBO joint
% Last edited 12/15/2021 by Lucien Peach

clear
close all
clc

addpath('DXFLib_v0.9.1')

% Specify inputs
r = 0.03; %[m]
n = 4;
nl = 4;
beta = pi/6; %[rad]
d0 = 0.05; 

% Outputs ls, l1, h0, dm, and alpha
[ls, l1, h0, alpha] = Origami_REBO_Parameters(r, n, beta, d0, nl);

% Create a figure that demonstrates the crease schematic:

% Though h1 and h2 will always be 0 for the purposes of these applications,
% they are included as inputs for the sake of consistency.
h1 = 0; %[m]
h2 = 0; %[m]

% Outputs graphing for REBO joint
[dataFoldREBO, m, lmax] = Origami_REBO_CreasePattern(n, nl, ls, l1, h1, h2, alpha);

% Create Duplication for Overlap Slide
[dataFoldNew] = StandaloneDuplication(dataFoldREBO, ls, n, lmax, 'REBO', h1);

% Plotting
figure()
hold on
for i = 1:size(dataFoldNew, 2)
    
    plot(dataFoldNew(i).x, dataFoldNew(i).y, 'color', dataFoldNew(i).color)    
    
end

% Label the plot for clarity
title({
    ('Standalone REBO Joint for Provided Parameters:')
    ['[r = ' num2str(r) ', n = ' num2str(n) ', beta = ' num2str(beta) ...
    ', nl = ' num2str(nl) ']']
    })

daspect([1 1 1])
axis off
set(gcf, 'color', 'w')
    
%     % Convert to DXF
    filename = (['REBO_r' num2str(r) '_n' num2str(n) '_l1_' num2str(l1) ...
        '.dxf']);
    GenerateDXF(filename, dataFoldREBO)
