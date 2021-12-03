% Creating a Crease Schematic - Origami prismatic joint
% Last edited 6/15/2021 by Lucien Peach
% 
% clear
% close all
% clc

% Specify inputs
r = 0.02; %[m]
n = 6;
nl = 6;
beta = pi/3; %[rad]
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
h1 = 0.05; %[m]
h2 = 0.05; %[m]

% h2 is defined as the region above the top orange line as in the Figure2.C
% depiction. See paper for reference.

% If no error is present (in requirement of l2), proceed with the
% graphing of the schematic and saving the file as a .dxf

% Outputs graphing for prismatic joint
[dataFoldE, m, lmax] = Origami_PrismaticJoint_CreasePattern(r, n, nl, ls, l1, dm, h0, h1, h2, alpha, beta);

% Plotting
figure()
hold on
for i = 1:size(dataFoldE, 2)
    
    plot(dataFoldE(i).x, dataFoldE(i).y, 'color', dataFoldE(i).color)    
    
end

% Label the plot for clarity
title({
    ('Origami Schematic C for Provided Parameters:')
    ['[r = ' num2str(r) ', n = ' num2str(n) ', beta = ' num2str(beta) ...
    ', h0 = ' num2str(h0) ', nl = ' num2str(nl) ', dm = ' num2str(l2) ']']
    })

daspect([1 1 1])
axis off
set(gcf, 'color', 'w')
    
%     % Convert to DXF
%     filename = (['FoldC_r' num2str(r) '_n' num2str(n) '_l1_' num2str(l1) ...
%         '_l2_' num2str(l2) '.dxf']);
%     GenerateDXF(filename, dataFoldE)
