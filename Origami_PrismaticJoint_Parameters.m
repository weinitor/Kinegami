function [ls, l1, h0, dm, alpha] = Origami_PrismaticJoint_Parameters(r, n, beta, d0, nl)
% ORIGAMI_PRISMATICJOINT_PARAMETERS - Finds the key variables needed to
% generate the crease pattern of the origami prismatic joint.

% Inputs:
%   r       - desired radius of folded origami linkage. 
%   n       - number of sides of folded origami linkage.
%   beta    - prismatic cone angle. 
%   d0      - zero configuration length.
%   nl      - number of folded prismatic layers.

% Outputs:
%   ls      - side length of folded origami linkage.
%   l1      - midsection height for single prismatic layer.
%   h0      - zero configuration length for a single prismatic layer.
%   dm      - the max length quantity, for use in determining height
%             boundaries of final origami schematic. 
%   alpha   - angular measurement for use in the generation of prismatic
%             fold region of schematic. 

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last edited 6/15/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


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

% Determine max length quantity
dm = nl * 2* l1;

% Display an error prompt if l2 is too small. Otherwise, return a value
% that can be checked to allow the code to proceed 
% if l2 > req
%     alert = 0;
% else
%     alert = msgbox('Choose a Larger Value for l2', 'Error', 'error');
% end

end
