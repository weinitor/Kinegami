function [r_mat] = RotationalMatrix(w, theta)
% ROTATIONALMATRIX - The rotational matrix about a unit vector w for theta
% degree.
%   w     - unit vector of a rotation axis.
%   theta - the angle of rotation (in rad).

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Last edited 7/27/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.txt for detail.


% The rotation axis must be a unit vector
w = w/norm(w);

wx = w(1);
wy = w(2);
wz = w(3);

c = cos(theta);
s = sin(theta);

% Initialize r_mat
r_mat = zeros(3, 3);

% Populate r_mat, the 3x3 rotational matrix

% Column 1
r_mat(1, 1) = (wx^2)*(1 - c) + c;
r_mat(2, 1) = wx*wy*(1 - c) + (wz*s);
r_mat(3, 1) = wx*wz*(1 - c) - (wy*s);

% Column 2
r_mat(1, 2) = wx*wy*(1 - c) - (wz*s);
r_mat(2, 2) = (wy^2)*(1 - c) + c;
r_mat(3, 2) = wy*wz*(1 - c) + (wx*s);

% Column 3
r_mat(1, 3) = wx*wz*(1 - c) + (wy*s);
r_mat(2, 3) = wy*wz*(1 - c) - (wx*s);
r_mat(3, 3) = (wz^2)*(1 - c) + c;

end
