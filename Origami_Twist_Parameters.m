function [x, l, ls] = Origami_Twist_Parameters(r, n, h, alpha)
% ORIGAMI_TWIST_PARAMETERS - Find the key variables needed to generate the
% crease pattern of the origami twist fitting.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last edited 6/22/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.txt for detail.


% Specify delta value for midsection height calculation
delta = pi*((n-2)/(2*n));

% Side length
ls = 2*r*cos(delta);

% Specify parameter for 1/2 base of midsection triangles
% x = r * mod(alpha, (2*pi)/n);
a = mod(alpha, (2*pi)/n);
x = ls/2 * (1-cos(a)+cot(pi/n)*sin(a));

% Length of midsection
l = sqrt((h^2) + (ls*csc(pi/n)*sin(pi/n-a/2)*sin(a/2))^2);

% this is how to calculate the marker
% m = floor(alpha/((2*pi)/n))*ls + x;

end
