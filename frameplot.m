function handle = frameplot(frame, color)
% FRAMEPLOT - Plot out the Op or Od frame.
% We know that the a vector is always normal to the frame, so we can plot
% the frame using this normal vector along with the centerpoint of the
% frame.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last edited 7/5/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


normal = frame(:, 1);
center = frame(:, 4);

for i = 1:3
    
    % Taking into account very small (but not 0) entries
    if normal(i) <= 10^-2 && normal(i) >= -10^-2
        
        normal(i) = 0;
        
    end
    
end

% Translating how to plot frames for each scenario
if normal(2) == 0 && normal(3) == 0 % Case for [1; 0; 0]
    
    % Create meshgrid with centerpoint adjustments
    [z, y] = meshgrid(-0.01:0.01:0.01);
    y = y + center(2);
    z = z + center(3);
    
    % x coordinates are all the center x coordinate
    x = zeros(size(y));
    x = x + center(1);
    
elseif normal(1) == 0 && normal(3) == 0 % Case for [0; 1; 0]
    
    % Create meshgrid with centerpoint adjustments
    [x, z] = meshgrid(-0.01:0.01:0.01);
    z = z + center(3);
    x = x + center(1);
    
    % y coordinates are all zero
    y = zeros(size(x));
    y = y + center(2);
    
elseif normal(1) == 0 && normal(2) == 0 % Case for [0; 0; 1]
    
    [x, y] = meshgrid(-0.01:0.01:0.01);
    x = x + center(1);
    y = y + center(2);
    
    % z coordinates are all zero
    z = zeros(size(x));
    z = z + center(3);
    
else % Case for any other scenario
    
    [x, y] = meshgrid(-0.01:0.01:0.01);
    
    A = normal(1);
    B = normal(2);
    C = normal(3);
    
    D = normal(1)*center(1) + normal(2)*center(2) + normal(3)*center(3);
    
    x = x + center(1);
    y = y + center(2);

    % Use plane formula to determine x values
    z = 1/C*(-A*x - B*y + D);  
          
end

% Plotting
handle = mesh(x, y, z, 'EdgeColor', 'none', 'FaceColor', color, 'FaceAlpha', 0.5);



end