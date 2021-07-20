% Self Assign
% Last Edited 7/20/2021 by Lucien Peach

function [TransformStruct] = SelfAssign(TransformStruct, r, n, JointStruct, N)

% Inputs: Necessary parameters as well as center locations and frame
% information for each "sphere"

% Assume value for beta to be constant
beta = pi/4; % [rad]

for i = 1:N+1
    
    % Assign r fields for sphere calculations
    
    % If revolute
    if JointStruct(i).type == 'R'
        
        TransformStruct(i).rs = r*sin(((n - 2)*pi) / (2*n))* ...
            tan(JointStruct(i).qm/ 4);
    
    % If prismatic
    elseif JointStruct(i).type == 'P'
        
        TransformStruct(i).rs = 1/4*JointStruct(i).q0*(2 + csc(beta));
    
    % Otherwise
    else
        
        TransformStruct(i).rs = r;
        
    end
    
end

for i = 1:N+1
    
    % Define the oinew field for later usage
    TransformStruct(i).oinew = TransformStruct(i).Oc(:, 4);    
    
end

red = [0.85, 0.325, 0.098];
yellow = [0.929, 0.694, 0.125];
green = [0.466, 0.674, 0.188];
blue = [0, 0.447, 0.741];

colorvector = [red; yellow; green; blue];

% Plot lines connecting consecutive spheres

% Initialize
xcenters = zeros(N+1, 1);
ycenters = zeros(N+1, 1);
zcenters = zeros(N+1, 1);

% Loop through to create centerpoint vectors
for i = 1:N+1
    
    xcenters(i) = TransformStruct(i).oinew(1);
    ycenters(i) = TransformStruct(i).oinew(2);
    zcenters(i) = TransformStruct(i).oinew(3);   
    
end

% We can now use the list of TransformStruct.Oc to provide us information
% on the new centroids of all of these spheres (joints)

% New visualization
figure()
set(gcf, 'color', 'w')
hold on

for i = 1:N+1  
    
    [TransformStruct(i).demonew] = SphericalSampling(TransformStruct(i).oinew, ...
    TransformStruct(i).rs, colorvector(i, :)); 
    
    % Plot vectors which demonstrate the new Oc 
    quiver3(TransformStruct(i).oinew(1), TransformStruct(i).oinew(2), ...
        TransformStruct(i).oinew(3), TransformStruct(i).Oc(1,1), ...
        TransformStruct(i).Oc(2,1), TransformStruct(i).Oc(3,1), ...
        'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'red');
    
    quiver3(TransformStruct(i).oinew(1), TransformStruct(i).oinew(2), ...
        TransformStruct(i).oinew(3), TransformStruct(i).Oc(1,2), ...
        TransformStruct(i).Oc(2,2), TransformStruct(i).Oc(3,2), ...
        'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'green');
    
    quiver3(TransformStruct(i).oinew(1), TransformStruct(i).oinew(2), ...
        TransformStruct(i).oinew(3), TransformStruct(i).Oc(1,3), ...
        TransformStruct(i).Oc(2,3), TransformStruct(i).Oc(3,3), ...
        'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'blue');
    
    grid on
    
end

% Plot
plot3(xcenters(:, 1), ycenters(:, 1), zcenters(:, 1), 'Color', 'k', ...
    'Linewidth', 4)

end