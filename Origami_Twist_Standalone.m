% Creating a Crease Schematic - Origami twist fitting
% Last edited 12/15/2021 by Lucien Peach

close all
clear 
clc

addpath('DXFLib_v0.9.1')

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
[dataFoldB, m, lmax] = Origami_Twist_CreasePattern(x, l, ls, n, h1, h2, r, h, alpha);

% Create Duplication for Overlap Slide
[dataFoldNew] = StandaloneDuplication(dataFoldB, ls, n, lmax);

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
% filename = (['Twist_r' num2str(r) '_n' num2str(n) '_h' num2str(h) ...
%     '_alpha' num2str(alpha) '.dxf']);
% GenerateDXF(filename, dataFoldB)
