% Creating a Crease Schematic - Fold Type B
% Last edited 6/7/2021 by Lucien Peach

clear
clc

% Specify inputs
r = 0.06; %[m]
n = 6;
h = 0.05; %[m]
alpha = 0; %[rad]

% Outputs values for x and l in [m]
[x, l, ls] = Origami_Twist_creasedesign(r, n, h, alpha);

% Create a figure that demonstrates the crease schematic

% Specify values for h1 and h2, the heights of the two tube sections
h1 = 0.03; %[m]
h2 = 0.03; %[m]

% Outputs midsection for graphing
[dataFoldB, m, lmax] = Origami_Twist_papercut(x, l, ls, n, h1, h2, r, h, alpha);
axis off

% Convert to DXF
% filename = (['FoldB_r' num2str(r) '_n' num2str(n) '_h' num2str(h) ...
%     '_alpha' num2str(alpha) '.dxf']);
% GenerateDXF(filename, dataFoldB)
