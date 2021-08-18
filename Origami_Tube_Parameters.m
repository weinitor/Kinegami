% Crease pattern parameters - Origami tube
% Last edited 6/9/2021 by Lucien Peach

% Function declaration
function [ls] = Origami_Tube_Parameters(r, n)

% Specify delta value for side length calculation and for use in parameter
% clarification
delta = pi*((n-2)/(2*n));

% Side length
ls = 2*r*cos(delta);

end
