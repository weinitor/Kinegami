% Basic Crease Design - Fold Type D (A.3) (Origami_RevoluteJoint)
% Last edited 6/11/2021 by Lucien Peach

% Function declaration
function [lengths, ls] = Origami_RevoluteJoint_creasedesign(r, n, theta_m)

% Specify delta value for side length and vertical lengths calculation
delta = pi*((n-2)/(2*n));

% Side length
ls = 2*r*cos(delta);

% Determine values of delta_i and l_i for each iteration. Store in arrays.
lengths = zeros(n, 1);
deltas = zeros(n, 1);

% Also create a vector for the sin^2(d_n) values as this will be populated
% iteratively
sinsquared_n = zeros(n, 1);
term2 = zeros(n, 1);


% Populate all values of delta_i, for use within l_i calculations
for i = 1:n
    deltas(i, 1) = pi * ((2 + n - (4*i)) / (2*n));
    
    % Use sin^2(delta_n) on each of the delta values for convenience
    sinsquared_n(i, 1) = sin(deltas(i, 1)) * sin(deltas(i, 1));
    
    % Determine denom of second term iteratively
%     term2(i, 1) = sin(pi/2 - deltas(i, 1) + delta) * sin(pi/2 - deltas(i, 1) + delta);
    
end

% Calculate values for sin^2(delta) and tan^2(theta_m / 4) prior to loop
sinsquaredelta = sin(delta) * sin(delta);
tansquare = tan(theta_m / 4) * tan(theta_m / 4);

% Determine lengths at each index for graphing
for j = 1:n
    lengths(j, 1) = r * sqrt((sinsquaredelta * tansquare) + ...
        sinsquared_n(j, 1));
end

end
