% Creating a Crease Schematic - Default origami tube
% Last edited 12/30/2021 by Lucien Peach

clear
clc

% Specify inputs
r = 0.02; %[m]
n = 4;

% Outputs ls, which will be used in graphing the default fold
[ls] = Origami_Tube_Parameters(r, n);

% Create a figure that demonstrates the crease schematic

% Specify value for h, the height of the tube section
h = 0.09; %[m]

% Outputs graphing for default tube
[dataFoldDefault, m, lmax] = Origami_Tube_CreasePattern(n, ls, h, r);

% Create Duplication for Overlap Slide
[dataFoldNew] = StandaloneDuplication(dataFoldDefault, ls, n, lmax, 'tube', h);

% Plotting
figure()
hold on
for i = 1:size(dataFoldNew, 2)
    
    plot(dataFoldNew(i).x, dataFoldNew(i).y, 'color', dataFoldNew(i).color)    
    
end

% Label the plot for clarity
title({
    ('Default Origami Schematic for Provided Parameters:')
    ['[r = ' num2str(r) ', n = ' num2str(n) ', h = ' num2str(h) ']']
    })

daspect([1 1 1])
axis off
set(gcf, 'color', 'w')

% Convert to DXF
filename = (['Tube_r' num2str(r) '_n' num2str(n) '_h' num2str(h) ...
    '.dxf']);
GenerateDXF(filename, dataFoldDefault)
