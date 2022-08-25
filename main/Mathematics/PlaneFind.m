function [a, b, c, d] = PlaneFind(plane)
% PLANEFIND - Simpler version of PlaneCheck that simply outputs values for
% a, b, c, and d such that the plane equation can be simply written.

% Inputs:
%   plane   - 3x2 matrix which includes information about the normal vector
%             and "center point" for the plane in question.

% Outputs:
%   a       - plane equation coefficient associated with x.
%   b       - plane equation coefficient associated with y.
%   c       - plane equation coefficient associated with z.
%   d       - constant contained in plane equation.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last edited 1/27/2022
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Plane represents the 3x2 matrix used to contain the normal vector and
% center point of the plane in question, while point is a point on the
% parallel z-axis we are seeking to analyze.
normal = plane(:, 1);
center = plane(:, 2);

% Determine coefficients for plane equation (a, b, c)
a = normal(1);
b = normal(2);
c = normal(3);

% Determine value of d by plugging info from center into plane equation
xa = center(1);
xb = center(2);
xc = center(3);

d = -(xa*a + xb*b + xc*c);

end