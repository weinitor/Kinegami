% Creating a Crease Schematic - Default origami tube
% Last edited 6/9/2021 by Lucien Peach

clear
clc

% Specify inputs
r = 0.2; %[m]
n = 6;

% Outputs ls, which will be used in graphing the default fold
[ls] = Origami_Tube_Parameters(r, n);

% Create a figure that demonstrates the crease schematic

% Specify value for h, the height of the tube section
h = 0.3; %[m]

% Outputs graphing for default tube
[dataFoldDefault, m, lmax] = Origami_Tube_CreasePattern(n, ls, h, r);
axis off

% Convert to DXF
% filename = (['FoldDefault_r' num2str(r) '_n' num2str(n) '_h' num2str(h) ...
%     '.dxf']);
% GenerateDXF(filename, dataFoldDefault)
