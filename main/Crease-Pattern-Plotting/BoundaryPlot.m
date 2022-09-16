function [dataFoldBoundary] = BoundaryPlot(n, ls, lmax_sum)
% BOUNDARYPLOT - Draws the boundary of the crease pattern.

% Inputs:
%   n                - number of sides of folded origami linkage.
%   ls               - side length of folded origami linkage. 
%   lmax_sum         - a cumulative counter that tracks the total height of
%                      all appended crease segments.

% Outputs:
%   dataFoldBoundary - data structure containing information relevant to
%                      the crease pattern boundary. 

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last Edited 6/29/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Counter used for data structure indexing
count = 1;

% Identify colors
black = [0, 0, 0];

% Boundary coordinates
boundary = [(n+1)*ls, 0; 0, 0; 0, lmax_sum; (n+1)*ls, lmax_sum; (n+1)*ls, 0];

dataFoldBoundary(count).x = boundary(:, 1);
dataFoldBoundary(count).y = boundary(:, 2);
dataFoldBoundary(count).color = black;

plot(dataFoldBoundary(count).x, dataFoldBoundary(count).y, 'color', ...
    dataFoldBoundary(count).color)
set(gcf, 'color', 'w')

daspect([1, 1, 1])
axis off

end