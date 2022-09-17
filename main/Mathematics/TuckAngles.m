function [tuckangle] = TuckAngles(r, n, phi, theta, split)
% TUCKANGLES - Calculates the angle that needs to be folded away for the
% hidden part of the elbow fitting. The inputs are the same as the elbow
% fitting, the output is an array of angles with the same size as
% "lengths". 

% Inputs:
%   r           - desired radius of folded origami linkage.
%   n           - number of sides of folded origami linkage.
%   phi         - rotational axis angle.
%   theta       - bending angle.
%   split       - determines whether theta splitting occurs.

% Outputs:
%   tuckangle   - an array of angles which can be used for generation of
%                 the origami elbow crease schematic.

% Authors: 
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Last edited 8/26/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


[lengths, ls] = Origami_Elbow_Parameters(r, n, phi, theta, split);
tuckangle = zeros(n, 1);
w = [cos(phi); sin(phi);0];
lengthsTemp = [lengths(n);lengths];
% the size of lengthsTemp is (n+2), its index is (n, 1, 2, ..., n, 1)

for k = 2:n+1
    % First find the angle around i vertex in the 2D crease pattern
    i = k;
    vfih = [-ls;lengthsTemp(i-1)-lengthsTemp(i);0];
    vfij = [ls;lengthsTemp(i+1)-lengthsTemp(i);0];
    thetafi = mod(SignedAngle(vfij, vfih, [0;0;1]),2*pi);
    
    % Then find the angle around i vertex on the truncated surface in 3D
    i = k-1;
    Ph = [r*cos(pi*(2*(i-1)-1)/n);r*sin(pi*(2*(i-1)-1)/n);0];
    Pi = [r*cos(pi*(2*i-1)/n);r*sin(pi*(2*i-1)/n);0];
    Pj = [r*cos(pi*(2*(i+1)-1)/n);r*sin(pi*(2*(i+1)-1)/n);0];
    v3dih = [Ph(1)-Pi(1);Ph(2)-Pi(2);vfih(2)];
    v3dij = [Pj(1)-Pi(1);Pj(2)-Pi(2);vfij(2)];
    an = RotationalMatrix(w,theta/2)*[0;0;1];
    theta3di = SignedAngle(v3dij, v3dih, an);
    
    tuckangle(i) = (thetafi-theta3di)/2;
end

% Outputs a size n+1 array of angles
tuckangle = [tuckangle;tuckangle(1)];

end