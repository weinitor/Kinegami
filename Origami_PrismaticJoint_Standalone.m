% ORIGAMI_PRISMATICJOINT_STANDALONE - Generates the crease pattern of the
% origami prismatic joint.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last Edited 12/15/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.txt for detail.


clear
close all
clc

% Add all the folders and subfolders to the search path
addpath(genpath(fileparts(mfilename('fullpath'))));

% Specify inputs
r = 0.03; %[m]
n = 4;
nl = 4;
beta = pi/6; %[rad]
d0 = 0.05; 
% h0 < rtan(beta)

% As a side note, provided the algorithms given, weird things happen if
% beta and n are both small. So not sure if some conditional should be
% implemented to prevent the program from running under these conditions

% Outputs ls, l1, alpha, and an alert message, if the parameter for l1 is
% not in accordance with algorithmic requirements
[ls, l1, h0, dm, alpha] = Origami_PrismaticJoint_Parameters(r, n, beta, d0, nl);

% req is output so that user can more easily adjust value of l2 to match
% given specifications for algorithm

% Create a figure that demonstrates the crease schematic

% Specify values for h1 and h2, the heights of the two tube sections
h1 = 0.02; %[m]
h2 = 0.02; %[m]

% h2 is defined as the region above the top orange line as in the Figure2.C
% depiction. See paper for reference.

% If no error is present (in requirement of l2), proceed with the
% graphing of the schematic and saving the file as a .dxf

% Outputs graphing for prismatic joint
[dataFoldE, m, lmax] = Origami_PrismaticJoint_CreasePattern(r, n, nl, ls, l1, dm, h0, h1, h2, alpha, beta);

% Create Duplication for Overlap Slide
[dataFoldNew] = StandaloneDuplication(dataFoldE, ls, n, lmax, 'prismatic', h1);

% Plotting
figure()
hold on
for i = 1:size(dataFoldNew, 2)
    
    plot(dataFoldNew(i).x, dataFoldNew(i).y, 'color', dataFoldNew(i).color)    
    
end

% Label the plot for clarity
title({
    ('Standalone Prismatic Joint for Provided Parameters:')
    ['[r = ' num2str(r) ', n = ' num2str(n) ', beta = ' num2str(beta) ...
    ', h0 = ' num2str(h0) ', nl = ' num2str(nl) ', dm = ' num2str(dm) ']']
    })

daspect([1 1 1])
axis off
set(gcf, 'color', 'w')
    
%     % Convert to DXF
    filename = (['Prismatic_r' num2str(r) '_n' num2str(n) '_l1_' num2str(l1) ...
        '_l2_' num2str(dm) '.dxf']);
    GenerateDXF(filename, dataFoldE)
