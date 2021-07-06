% Joint Assignment
% Last Edited 6/25/2021 by Lucien Peach

% Make r = 0.01
% Make q0 = {0; 0; 0; 0.05; 0}
% Make qm = {pi; pi; pi; 0; 0}
% Both of the above will be fields of JointStruct
% Make JointStruct.type = {R; R; R; P; 0}

function [TransformStruct] = JointAssignment(D, r, n, JointStruct, N, theta_mod, fingertip)

TransformStruct(N+1) = struct();

for i = (N+1):-1:1
    
    % The i index refers to the lower right value for regular and upper
    % left value for inverse
    TransformStruct(i).T = HomogeneousTransform(i, D);
    TransformStruct(i).inverse = InverseHomogeneousTransform(i, D);
    
    % This will be used to find the value of 0{T}N+1
    if i == (N+1)
        TransformStruct(N+1).net = TransformStruct(i).T;
    else
        TransformStruct(N+1).net = TransformStruct(i).T * TransformStruct(N+1).net;
    end
    
end

% Extract O(N+1) and Oc(N+1) from NetTransform (4x4 matrix)
Ox = TransformStruct(N+1).net(1:3, 1);
Oy = TransformStruct(N+1).net(1:3, 2);
Oz = TransformStruct(N+1).net(1:3, 3);
Oc = TransformStruct(N+1).net(1:3, 4);

% O and Oc, initial axis of translation
TransformStruct(N+1).O = [Ox, Oy, Oz, Oc];

if strcmp(fingertip, 'z')
    TransformStruct(N+1).Oc = [Oz, Ox, Oy, Oc];
elseif strcmp(fingertip, 'x')
    TransformStruct(N+1).Oc = [Ox, Oy, Oz, Oc];
elseif strcmp(fingertip, 'y')
    TransformStruct(N+1).Oc = [Oy, Oz, Ox, Oc];
end

TransformStruct(N+1).zaxis = Oz.';
TransformStruct(N+1).xaxis = Ox.';

% Define r and oi for initial sphere (recall .oi etc. is row)
TransformStruct(N+1).oi = Oc.';
TransformStruct(N+1).rs = r;

% Assume value for beta to be constant
beta = pi/4; % [rad]

for i = 1:N
    
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

% Loop through homogeneous transforms
for i = N:-1:1
    
    % Multiply inverse matrix and standard homogeneous to compute net
    TransformStruct(i).net = TransformStruct(i+1).net * ...
        TransformStruct(i+1).inverse;
    
    % Extract Oi from TransformStruct(i).net
    Ox = TransformStruct(i).net(1:3, 1);
    Oy = TransformStruct(i).net(1:3, 2);
    Oz = TransformStruct(i).net(1:3, 3);
    Oc = TransformStruct(i).net(1:3, 4);

    % O Computation
    TransformStruct(i).O = [Ox, Oy, Oz, Oc];
    
    % Translation Axis
    TransformStruct(i).zaxis = Oz.';
    
    % X Axis
    TransformStruct(i).xaxis = Ox.';
    
    % Define r and oi for spheres
    TransformStruct(i).oi = Oc.';

end

% Initial visualization of spheres
figure()
set(gcf, 'color', 'w')
hold on

for i = 1:N+1
    
    [TransformStruct(i).demo] = SphericalSampling(TransformStruct(i).oi, ...
        TransformStruct(i).rs, 'none'); 
end

% Sphere bounding fundamentals for (N+1) ---------------------------

% Sphere bounding data
figure()
set(gcf, 'color', 'w')
hold on

[TransformStruct(N+1).demo] = SphericalSampling(TransformStruct(N+1).oi, ...
    TransformStruct(N+1).rs, 'none');

% Concatenate to express array of spheres
TransformStruct(N+1).democumul = TransformStruct(N+1).demo;

% Find new sphere center and store to index-1.optimized
[R, C, Xb] = ExactMinBoundSphere3D(TransformStruct(N+1).democumul);

TransformStruct(N+1).Xb = Xb;
TransformStruct(N+1).oinew = C;
TransformStruct(N+1).rnew = R;

% Output new plot of cumulative sphere
[TransformStruct(N+1).democumul] = SphericalSampling(TransformStruct(N+1).oinew, ...
    TransformStruct(N+1).rnew, 'none');

TransformStruct(N+1).oilarge = TransformStruct(N+1).oi;

% Loop through sphere bounding and minimization (Fun!)
for i = (N+1):-1:2
    
    [TransformStruct] = SphereMinimization(TransformStruct, i, r, 'none');  
     
end

% Final sphere ("joint") visualization
figure()
set(gcf, 'color', 'w')
hold on

red = [0.85, 0.325, 0.098];
yellow = [0.929, 0.694, 0.125];
green = [0.466, 0.674, 0.188];
blue = [0, 0.447, 0.741];

colorvector = [red; yellow; green; blue];

for i = 1:N+1
    
    % Original spheres
%     [TransformStruct(i).demo] = SphericalSampling(TransformStruct(i).oi, ...
%     TransformStruct(i).r, 0.7*colorvector(i, :)); 
% 
%     quiver3(TransformStruct(i).oi(1), TransformStruct(i).oi(2), ...
%         TransformStruct(i).oi(3), TransformStruct(i).axis(1), ...
%         TransformStruct(i).axis(2), TransformStruct(i).axis(3), ...
%         'AutoScaleFactor', 1, 'Linewidth', 1)
    
    % Plot final new sphere locations
    [TransformStruct(i).demonew] = SphericalSampling(TransformStruct(i).oinew, ...
        TransformStruct(i).rs, colorvector(i, :)); 
    
    
    % Plot vectors which demonstrate the axis along which each sphere is
    % restricted to move
    quiver3(TransformStruct(i).oinew(1), TransformStruct(i).oinew(2), ...
        TransformStruct(i).oinew(3), TransformStruct(i).zaxis(1), ...
        TransformStruct(i).zaxis(2), TransformStruct(i).zaxis(3), ...
        'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 0.8*colorvector(i, :));
    
    quiver3(TransformStruct(i).oinew(1), TransformStruct(i).oinew(2), ...
        TransformStruct(i).oinew(3), TransformStruct(i).xaxis(1), ...
        TransformStruct(i).xaxis(2), TransformStruct(i).xaxis(3), ...
        'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'k');
    
    grid on
end

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

% Plot
plot3(xcenters(:, 1), ycenters(:, 1), zcenters(:, 1), 'Color', 'k', ...
    'Linewidth', 4)
    

% Loop through joint coordinate reassignment
for i = 1:N
    
    % Redefine x y and z, as well as centroid, based on sphere reassignment
    Ox = TransformStruct(i).net(1:3, 1);
    Oy = TransformStruct(i).net(1:3, 2);
    Oz = TransformStruct(i).net(1:3, 3);
    Oc = TransformStruct(i).oinew.';
    
    % If revolute
    if JointStruct(i).type == 'R'
        
        ai = rotationalmatrix(Oz, theta_mod(i))*Ox;
        bi = rotationalmatrix(Oz, theta_mod(i))*Oy;
        
        TransformStruct(i).Oc = [ai, bi, Oz, Oc];
    
    % If prismatic
    elseif JointStruct(i).type == 'P'
        
        TransformStruct(i).Oc = [Oz, Ox, Oy, Oc];      
    
    % Other joints
    else
        
        TransformStruct(i).Oc = [Ox, Oy, Oz, Oc];
        
    end
    
end

% We can now use the list of TransformStruct.Oc to provide us information
% on the new centroids of all of these spheres (joints)

% New visualization
hold off
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