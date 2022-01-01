% Self Assign
% Last Edited 7/20/2021 by Lucien Peach

function [TransformStruct] = SelfAssign(TransformStruct, r, n, JointStruct, N, plotoption)

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
    
    % Otherwise, for fingertip
    elseif JointStruct(i).type == 'F'
        
        TransformStruct(i).rs = r*sin(((n - 2)*pi) / (2*n))* ...
            tan(JointStruct(i).qm/ 4);
     
    % For V Joint
    else 
        
        TransformStruct(i).rs = 0;
           
    end
    
end

for i = 1:N+1
    
    % Define the oinew field for later usage
    TransformStruct(i).oinew = TransformStruct(i).Oc(:, 4);  
    
    if i == N+1
        
        % Define oilarge for initial sphere (fingertip)
        TransformStruct(i).oilarge = TransformStruct(i).oinew;        
        
    end
    
end

if strcmp(plotoption, 'on') == 1
    figure()
    set(gcf, 'color', 'w')
    hold on
end

% Add section to enable individual sphere and cumulative sphere plotting
for i = (N+1):-1:2
    
    for j = size(TransformStruct, 2):-1:i-1

        [TransformStruct(j).demo] = SphericalSampling(TransformStruct(j).oinew, ...
            TransformStruct(j).rs, 'none', plotoption);

        % If on the last index, the cumulative sphere points is just the final
        % sphere itself
        if j == size(TransformStruct, 2)

            TransformStruct(j).democumul = TransformStruct(j).demo;  
        else

            TransformStruct(j).democumul = [TransformStruct(j).demo; TransformStruct(j+1).democumul];
        end
    
    end
    
    if i == N+1
        
        % Define overlap sphere for fingertip
        TransformStruct(i).rnew = TransformStruct(i).rs;
        
        [TransformStruct(i).spheresnet] = SphericalSampling(TransformStruct(i).oilarge, ...
            TransformStruct(i).rnew, 'none', plotoption);
         
    end

    % Find new sphere center and store to index-1.optimized
    [R, C, Xb] = ExactMinBoundSphere3D(TransformStruct(i-1).democumul);

    TransformStruct(i-1).Xb = Xb;
    TransformStruct(i-1).oilarge = C;
    TransformStruct(i-1).rnew = R;

    % Output new plot of cumulative sphere
    [TransformStruct(i-1).spheresnet] = SphericalSampling(TransformStruct(i-1).oilarge, ...
        TransformStruct(i-1).rnew, 'none', plotoption);  
    
end

if strcmp(plotoption, 'on') == 1
    hold off
end

% After Intersection runs, oinew locations will be indicative of position
% with no path overlaps

% Sphere color assignment
black = [0,0,0];
red = [0.6350, 0.0780, 0.1840];
orange = [0.8500, 0.3250, 0.0980];
yellow = [0.929, 0.694, 0.125];
green = [0.466, 0.674, 0.188];
blue = [0, 0.447, 0.741];
purple = [0.4940, 0.1840, 0.5560];

colorvector = [black; red; orange; yellow; green; blue; purple];

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
if strcmp(plotoption, 'on') == 1
    figure()
    set(gcf, 'color', 'w')
    hold on
end

for i = 1:N+1  
    
    [TransformStruct(i).demonew] = SphericalSampling(TransformStruct(i).oinew, ...
    TransformStruct(i).rs, colorvector(i, :), plotoption); 
    
    if strcmp(plotoption, 'on') == 1
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
    
end

% Plot
if strcmp(plotoption, 'on') == 1 
    plot3(xcenters(:, 1), ycenters(:, 1), zcenters(:, 1), 'Color', 'k', ...
        'Linewidth', 4)
end

end