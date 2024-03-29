function [lengths, ls] = Origami_RevoluteJoint_Parameters(r, n, theta_m)
% ORIGAMI_REVOLUTEJOINT_PARAMETERS - Find the key variables needed to
% generate the crease pattern of the origami revolute joint.

% Inputs:
%   r           - desired radius of folded origami linkage.
%   n           - number of sides of folded origami linkage.
%   theta_m     - total bending angle of revolute joint.

% Outputs:
%   lengths     - vector of length measurements for use in schematic
%                 plotting.
%   ls          - side length of folded origami linkage.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last edited 6/17/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Specify delta value for side length and vertical lengths calculation
delta = pi*((n-2)/(2*n));

% Side length
ls = 2*r*cos(delta);

% Find the half joint height h, column height or max(li) lmax, alpha
h = r*sin(delta)*tan(theta_m/4);
lmax = r*sin(delta)/cos(theta_m/4);
alpha = atan(ls/(2*lmax));

% Determine values of delta_i and l_i for each iteration. Store in arrays.
lengths = zeros(n, 1);
deltas = zeros(n, 1);
gammas_n = zeros(n, 1);

% Also create a vector for the sin^2(d_n) values as this will be populated
% iteratively
lfront_n = zeros(n, 1);
den_n = zeros(n, 1);


% Populate all values of delta_i, for use within l_i calculations
for i = 1:n
    deltas(i, 1) = pi * ((2 + n - (4*i)) / (2*n));
    
    % lfront
    lfront_n(i, 1) = sqrt(h^2 + r^2*sin(deltas(i, 1))*sin(deltas(i, 1)));
    
    % Gamma
    gammas_n(i, 1) = atan(lfront_n(i, 1)/(r*abs(cos(deltas(i, 1)))));
    
    % Denominator
    den_n(i, 1) = cos(pi/2 -gammas_n(i, 1) - alpha);
    
end

% Determine lengths at each index for graphing
for j = 1:n
    lengths(j, 1) = lfront_n(j, 1) / den_n(j, 1);
end

end
