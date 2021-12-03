% Creating a Crease Schematic - Origami twist fitting
% Last edited 6/7/2021 by Lucien Peach

clear
clc

% Specify inputs
r = 0.06; %[m]
n = 6;
h = 0.05; %[m]
alpha = pi/4; %[rad]

% Outputs values for x and l in [m]
[x, l, ls] = Origami_Twist_Parameters(r, n, h, alpha);

% Create a figure that demonstrates the crease schematic

% Specify values for h1 and h2, the heights of the two tube sections
h1 = 0.03; %[m]
h2 = 0.03; %[m]

% Outputs midsection for graphing
[dataFoldB, m, lmax] = Origami_Twist_CreasePattern(x, l, ls, n, h1, h2, r, h, alpha);

% Plotting
figure()
hold on
for i = 1:size(dataFoldD, 2)
    
    plot(dataFoldB(i).x, dataFoldB(i).y, 'color', dataFoldB(i).color)    
    
end

% Label the plot for clarity
title({
    ('Origami Schematic B for Provided Parameters:')
    ['[r = ' num2str(r) ', n = ' num2str(n) ', h = ' num2str(h) ', alpha = ' num2str(alpha) ']']
    })

daspect([1 1 1])
axis off
set(gcf, 'color', 'w')

% Convert to DXF
% filename = (['FoldB_r' num2str(r) '_n' num2str(n) '_h' num2str(h) ...
%     '_alpha' num2str(alpha) '.dxf']);
% GenerateDXF(filename, dataFoldB)
