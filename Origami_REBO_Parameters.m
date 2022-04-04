% Crease pattern parameters - Origami REBO joint
% Last edited 4/3/2022 by Lucien Peach

% Function declaration
function [ls, l1, h0, alpha] = Origami_REBO_Parameters(r, n, beta, d0, nl)

% Specify delta value for side length calculation
delta = pi*((n-2)/(2*n));

% Side length
ls = 2*r*cos(delta);

% Determine value of alpha for use in angular measurements
alpha = (pi - (2*pi*cos(beta)/n))/2;

% Determine h0
h0 = d0/nl;

% Determine value of l1, the height of the midsection 
l1 = h0*csc(beta)/2;

end
