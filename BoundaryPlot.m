function [dataFoldBoundary] = BoundaryPlot(n, ls, lmax_sum)
% BOUNDARYPLOT - Draw the boundary of the crease pattern.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last Edited 6/29/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.txt for detail.


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