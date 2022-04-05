function [ls] = Origami_Tube_Parameters(r, n)
% ORIGAMI_TUBE_PARAMETERS - Find the key variables needed to generate the
% crease pattern of the origami prism tube.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last edited 6/9/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Specify delta value for side length calculation and for use in parameter
% clarification
delta = pi*((n-2)/(2*n));

% Side length
ls = 2*r*cos(delta);

end
