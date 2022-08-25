function [Concat, handle] = SphericalSampling(oi, r, colorinput, plotoption)
% SPHERICALSAMPLING - Provides visualization of spherical region based on
% center point, radius of sphere, and other minor settings.

% Inputs:
%   oi          - the component of Oi which describes the location of the
%                 point in space.
%   r           - the radius unadjusted by any joint parameters (pure r).
%   colorinput  - string input which determines the color of the plotted
%                 sphere, if so desired.
%   plotoption  - string input which dictates plotting.

% Outputs:
%   Concat      - matrix of values which contains 3 columns of sphere data
%                 for use within ExactMinBoundSphere3D.m
%   handle      - plot handle which allows for future edits to plot
%                 settings or to call the plot if needed.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Created 6/21/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Specify colorinput as 'none' if random color is desired

% oi is the component of Oi that describes the location of the point in
% space. For construction of spherical units, axes orientation is not
% necessary (will be necessary for comparison and translation of spheres)

% r is the common radius of all joint components

% Create values for azimuth (TH) and elevation (PHI)
TH = 2*pi*rand(1, 1e4);
PHI = asin(-1 + 2*rand(1, 1e4));
R = r;

% Extract X Y and Z 
[X, Y, Z] = sph2cart(TH, PHI, R);

% Convert row vectors to column vectors. Add offset for origin location.
% Column vector format will be needed for ExactMinBoundSphere3D
X = X.' + oi(1);
Y = Y.' + oi(2);
Z = Z.' + oi(3);

% Concatenate for ExactMinBoundSphere3D
Concat = [X, Y, Z];

TF = strcmp(colorinput, 'none');

if TF == 1
    
    if strcmp(plotoption, 'on') == 1
        % 3D Plot (can comment out if no visualization needed)
        handle = plot3(X, Y, Z, '.', 'markersize', 1);
        axis equal vis3d

        xlabel('x')
        ylabel('y')
        zlabel('z')
    end

else
    
    if strcmp(plotoption, 'on') == 1
        % 3D Plot (can comment out if no visualization needed)
        handle = plot3(X, Y, Z, '.', 'markersize', 1, 'color', colorinput);
        axis equal vis3d

        xlabel('x')
        ylabel('y')
        zlabel('z')
    end
 
end


end