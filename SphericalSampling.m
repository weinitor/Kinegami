% Spherical Sampling
% Last Edited 6/23/2021 by Lucien Peach

function [Concat] = SphericalSampling(oi, r, colorinput)

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
    
    % 3D Plot (can comment out if no visualization needed)
    plot3(X, Y, Z, '.', 'markersize', 1)
    axis equal vis3d
     
    xlabel('x')
    ylabel('y')
    zlabel('z')

else
   % 3D Plot (can comment out if no visualization needed)
    plot3(X, Y, Z, '.', 'markersize', 1, 'color', colorinput)
    axis equal vis3d
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
 
end


end