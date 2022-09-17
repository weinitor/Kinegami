function [theta] = SignedAngle(a, b, n)
% SIGNEDANGLE - Finding the signed angle between two vectors "a" and "b"
% with respect to vector "a", around the normal vector "n". Outputs an
% angle "theta".

% Inputs:
%   a       - vector "a". See function description for details.
%   b       - vector "b". See function description for details.
%   n       - normal vector for utilization in cross product.

% Outputs:
%   theta   - output signed angle, in radians.

% Authors: 
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Last edited 8/6/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Make sure all vectors are normalised
a = a/norm(a);
b = b/norm(b);
n = n/norm(n);

% find the signed angle
theta = atan2(dot(cross(a,b),n),dot(a,b));

end

