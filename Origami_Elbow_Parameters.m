% Crease pattern parameters - Origami elbow fitting
% Last edited 8/14/2021 by Lucien Peach

% Function declaration
function [lengths, ls] = Origami_Elbow_Parameters(r, n, phi, theta, split)

% Flip the direction is theta is less tha 0
if theta < 0
    
    theta = abs(theta);
    phi = pi + phi;
    
end

if theta > pi/2 && strcmp(split, 'off') ~= 1
    
    theta = theta/2;
    
end

% If theta is sufficiently small
if theta <= 0.015
    
    lengths = zeros(n+1, 1);
   
    % Specify delta value for side length calculation
    delta = pi*((n-2)/(2*n));

    % Side length
    ls = 2*r*cos(delta);
    
end

% Otherwise
if theta > 0.015

    % Determine values of delta_i and l_i for each iteration. Store in arrays.

    lengths = zeros(n+1, 1);
    deltas = zeros(n, 1);
    sindeltas = zeros(n, 1);

    % Populate all values of delta_i so that elbow fitting offset can be found
    % for use in populating the lengths array
    for i = 1:n
        deltas(i, 1) = - phi + (pi * ((2*i - 1) / n));

        % Use sin operation on each of the delta values
        sindeltas(i, 1) = sin(deltas(i, 1))*tan(theta/2);    

    end

    % Use this minimum value to find the elbow fitting offset
    d = r * tan(theta/2);

    % Determine lengths at each index for graphing
    for j = 1:n
        lengths(j, 1) = (r * (sindeltas(j, 1))) + d;
    end

    % Duplicate value of initial value and append to end for glue region
    lengths(end, 1) = lengths(1, 1);

    % Specify delta value for side length calculation
    delta = pi*((n-2)/(2*n));

    % Side length
    ls = 2*r*cos(delta);
    
end

end
