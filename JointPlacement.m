% Joint Placement
% Last Edited 1/28/2021 by Lucien Peach

% Make r = 0.01
% Make q0 = {0; 0; 0; 0.05; 0}
% Make qm = {pi; pi; pi; 0; 0}
% Both of the above will be fields of JointStruct
% Make JointStruct.type = {R; R; R; P; 0}

function [TransformStruct] = JointPlacement(D, r, n, JointStruct, N, theta_mod, fingertip, plotoption)

% Does this process remain for JointPlacement? Unsure but kept for now.
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

for i = 1:N+1
    
    % Assign r fields for sphere calculations
    
    % If revolute
    if JointStruct(i).type == 'R'
        
        TransformStruct(i).rs = r*sin(((n - 2)*pi) / (2*n))* ...
            tan(JointStruct(i).qm/ 4);
    
    % If prismatic
    elseif JointStruct(i).type == 'P'
        
        TransformStruct(i).rs = 1/4*JointStruct(i).q0*(2 + csc(beta));
    
    % If waypoint
    elseif JointStruct(i).type == 'W'
        
        TransformStruct(i).rs = 0;
    
    % Otherwise, for fingetip
    else
        
        TransformStruct(i).rs = r*sin(((n - 2)*pi) / (2*n))* ...
            tan(JointStruct(i).qm/4);
        
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
if strcmp(plotoption, 'on') == 1
    figure()
    set(gcf, 'color', 'w')
    hold on
    grid on
end

% Initial Sphere Visualization (will plot if plotoption is true)
for i = 1:N+1
    
    [TransformStruct(i).oidata] = SphericalSampling(TransformStruct(i).oi, ...
        TransformStruct(i).rs, 'none', plotoption);
    
end

% Initial bounding sphere is same as first sphere
TransformStruct(N+1).oib = TransformStruct(N+1).oi;
TransformStruct(N+1).rb = TransformStruct(N+1).rs;
TransformStruct(N+1).boundingdata = TransformStruct(N+1).oidata;

% Assignment Looping
for i = N+1:-1:2
    
    % Normal vector adjustment (column vector)
    TransformStruct(i).normal = -TransformStruct(i).Oc(:, 1);

    % Plane generation for tangent plane
    TransformStruct(i).tangent(:, 4) = TransformStruct(i).oib.' + ...
      TransformStruct(i).rb * TransformStruct(i).normal;

    TransformStruct(i).tangent(:, 1) = TransformStruct(i).normal;

    % Storing plane data
    TransformStruct(i).P1 = zeros(3, 2);
    TransformStruct(i).P1(:, 1) = TransformStruct(i).tangent(:, 1);
    TransformStruct(i).P1(:, 2) = TransformStruct(i).tangent(:, 4);  

    % Plane generation for parallel plane
    TransformStruct(i).parallel(:, 4) = TransformStruct(i).oib.' + ...
      (TransformStruct(i).rb + TransformStruct(i-1).rs + 4*r) * ...
      TransformStruct(i).normal;

    TransformStruct(i).parallel(:, 1) = TransformStruct(i).normal;

    % Storing parallel plane data
    TransformStruct(i).P1par = zeros(3, 2);
    TransformStruct(i).P1par(:, 1) = TransformStruct(i).parallel(:, 1);
    TransformStruct(i).P1par(:, 2) = TransformStruct(i).parallel(:, 4);  

    % Determine location of first waypoint
    [waypoint1] = IntersectionSolver(TransformStruct(i).P1, ...
      TransformStruct(i).oi.', TransformStruct(i).normal);

    % Store to structure as waypoint1
    TransformStruct(i).waypoint1 = [TransformStruct(i).Oc(:, 1), ...
      TransformStruct(i).Oc(:, 2), TransformStruct(i).Oc(:, 3), ...
      waypoint1];

    % Determine intersection of preceding z-axis (for sphere being placed)
    % and plane, if it occurs.

    % Run PlaneCheck, which will determine if the z-axis is on the far or
    % near side of the parallel plane.
    if dot(TransformStruct(i).normal, TransformStruct(i-1).zaxis) == 0

      planeval = PlaneCheck(TransformStruct(i).P1par, TransformStruct(i-1).oi);   

    end

    % For parallel, conflicting case
    if dot(TransformStruct(i).normal, TransformStruct(i-1).zaxis) == 0 && ...
            planeval < 0

        % Assign normal vector
        TransformStruct(i).normal2 = -TransformStruct(i-1).zaxis;

        % Assign new planes
        TransformStruct(i).tangent2(:, 4) = TransformStruct(i).oib.' + ...
            TransformStruct(i).rb * TransformStruct(i).normal2;

        TransformStruct(i).tangent2(:, 1) = TransformStruct(i).normal2;

        % Storing plane data
        TransformStruct(i).P2 = zeros(3, 2);
        TransformStruct(i).P2(:, 1) = TransformStruct(i).tangent2(:, 1);
        TransformStruct(i).P2(:, 2) = TransformStruct(i).tangent2(:, 4);  

        % Plane generation for parallel plane
        TransformStruct(i).parallel2(:, 4) = TransformStruct(i).oib.' + ...
            (TransformStruct(i).rb + TransformStruct(i-1).rs + 4*r) * ...
            TransformStruct(i).normal2;

        TransformStruct(i).parallel2(:, 1) = TransformStruct(i).normal2;

        % Storing parallel plane data in useful format
        TransformStruct(i).P2par = zeros(3, 2);
        TransformStruct(i).P2par(:, 1) = TransformStruct(i).parallel2(:, 1);
        TransformStruct(i).P2par(:, 2) = TransformStruct(i).parallel2(:, 4);  

        % Determining waypoint2
        ref2 = TransformStruct(i).waypoint1 + r*TransformStruct(i).normal;
        [waypoint2] = IntersectionSolver(TransformStruct(i).P2, ...
            ref2, TransformStruct(i).normal2);

        % Rotational Frame Transform
        rotationvec = cross(TransformStruct(i).normal, TransformStruct(i).normal2);
        rotationvec = rotationvec/norm(rotationvec);
      
        [rotation_matrix] = RotationalMatrix(rotationvec, pi/2);
      
        OcMatrix = [TransformStruct(i).Oc(:, 1), TransformStruct(i).Oc(:, 2), ...
            TransformStruct(i).Oc(:, 3)];
        ORotation = OcMatrix * rotation_matrix;
      
        % Store to structure as waypoint2
        TransformStruct(i).waypoint2 = [ORotation(:, 1), ORotation(:, 2), ...
            ORotation(:, 3), waypoint2];

    % For non-parallel case, non-conflicting case
    else
        
        % Assign n2 vector
        TransformStruct(i).normal2 = TransformStruct(i).normal;

        % Assign new planes
        TransformStruct(i).P2 = TransformStruct(i).P1;
        TransformStruct(i).P2par = TransformStruct(i).P1par;

        % Assign waypoint2
        TransformStruct(i).waypoint2 = TransformStruct(i).waypoint1;

    end

    % Sphere bounding process

    % Oc assignment
    % Joint Type Dictates Ox, Oy, Oz, Oc order. Clearly define these.
    Ox = TransformStruct(i-1).O(:, 1);
    Oy = TransformStruct(i-1).O(:, 2);
    Oz = TransformStruct(i-1).O(:, 3);
    Oc = TransformStruct(i-1).O(:, 4);
    
    % Since fingertip case (N+1) is already addressed in lines 29-44,
    % worry only about Prismatic and Revolute cases
    if JointStruct(i-1).type == 'R'
        
        TransformStruct(i-1).Oc = [Ox, Oy, Oz, Oc];
        
    elseif JointStruct(i-1).type == 'P'
        
        TransformStruct(i-1).Oc = [Oz, Ox, Oy, Oc];        
                
    end
    
    % Now that we have bare bones template for Oc, vary values
    % Begin by determining u and t    
    if -TransformStruct(i-1).Oc(:, 1).' * TransformStruct(i).normal2 >= 0
        
        u = 1;       
        
    elseif -TransformStruct(i-1).Oc(:, 1).' * TransformStruct(i).normal2 < 0
        
        u = -1;
        
    end
    
    % Optimization based on constraints to determine value of t
    
    % Statement for minimization
    obj = @(delta) norm(TransformStruct(i-1).Oc(:, 4) + delta*Oz - ...
        TransformStruct(i).waypoint2(:, 4));

    % Initial Guess
    delta0 = 1;
    
    % Run PlaneFind to determine variables for a, b, c, d for plane in
    % question.    
    [a, b, c, d] = PlaneFind(TransformStruct(i).P2par);

    % Inequality requirement:
    % a*oi(1) + b*oi(2) + c*oi(3) + d >= 0

    % Express other parameters as empty cells (not used)
    A = [];
    B = [];
    Aeq = [];
    Beq = [];
    LB = [];
    UB = [];
    nonlincon = @nlcon;

    % Output optimal value for delta (or t)
    t_val = fmincon(obj, delta0, A, B, Aeq, Beq, LB, UB, nonlincon);
    
    % Now, use values of u and t to find Oc(i-1)
    if JointStruct(i-1).type == 'R'
        
        TransformStruct(i-1).Oc = [u*Ox, u*Oy, Oz, Oc + t_val*Oz];
        
    elseif JointStruct(i-1).type == 'P'
        
        TransformStruct(i-1).Oc = [u*Oz, u*Ox, Oy, Oc + t_val*Oz];        
                
    end
    
    % Create O DataSet
    % Line 22 of Script (why necessary?)
    
    % New sphere creation
    [TransformStruct(i-1).oidata] = SphericalSampling(TransformStruct(i-1).Oc(:, 4), ...
        TransformStruct(i-1).rs, 'none', plotoption);
    
    % Bounding sphere data
    TransformStruct(i-1).boundingdata = [TransformStruct(i).boundingdata; ...
        TransformStruct(i-1).oidata];
    
    % Find new oib and rb values
    [R, C, Xb] = ExactMinBoundSphere3D(TransformStruct(i-1).boundingdata);
    
    TransformStruct(i-1).oib = C;
    TransformStruct(i-1).rb = R;
    TransformStruct(i-1).Xb = Xb;
    
    % Output plot of new bounding sphere
    [TransformStruct(i-1).boundingsummary] = SphericalSampling(TransformStruct(i-1).oib, ...
        TransformStruct(i-1).rb, 'none', plotoption);
    
end

% Nonlinear Constraint Function Call
function [constraint, ceq] = nlcon(delta)
    constraint = 0 - (a*(TransformStruct(i-1).Oc(1, 4) + delta*Oz(1)) + ...
        b*(TransformStruct(i-1).Oc(2, 4) + delta*Oz(2)) + ...
        c*(TransformStruct(i-1).Oc(3, 4) + delta*Oz(3)) + d);
    ceq = [];
end

% Final sphere ("joint") visualization
if strcmp(plotoption, 'on') == 1
    figure()
    set(gcf, 'color', 'w')
    hold on
end

black = [0,0,0];
red = [0.6350, 0.0780, 0.1840];
orange = [0.8500, 0.3250, 0.0980];
yellow = [0.929, 0.694, 0.125];
green = [0.466, 0.674, 0.188];
blue = [0, 0.447, 0.741];
purple = [0.4940, 0.1840, 0.5560];

colorvector = [black; red; orange; yellow; green; blue; purple];

for i = 1:N+1
        
    % Plot final new sphere locations
    [TransformStruct(i).demonew] = SphericalSampling(TransformStruct(i).Oc(:, 4), ...
        TransformStruct(i).rs, colorvector(i, :), plotoption); 
       
    % Plot vectors which demonstrate the axis along which each sphere is
    % restricted to move
    if strcmp(plotoption, 'on') == 1
        quiver3(TransformStruct(i).Oc(1, 4), TransformStruct(i).Oc(2, 4), ...
            TransformStruct(i).Oc(3, 4), TransformStruct(i).zaxis(1), ...
            TransformStruct(i).zaxis(2), TransformStruct(i).zaxis(3), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 0.8*colorvector(i, :));

        quiver3(TransformStruct(i).Oc(1, 4), TransformStruct(i).Oc(2, 4), ...
            TransformStruct(i).Oc(3, 4), TransformStruct(i).xaxis(1), ...
            TransformStruct(i).xaxis(2), TransformStruct(i).xaxis(3), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'k');

        grid on
    end
end

% Plot lines connecting consecutive spheres

% Initialize
xcenters = zeros(N+1, 1);
ycenters = zeros(N+1, 1);
zcenters = zeros(N+1, 1);

% Loop through to create centerpoint vectors
for i = 1:N+1
    
    xcenters(i) = TransformStruct(i).Oc(1, 4);
    ycenters(i) = TransformStruct(i).Oc(2, 4);
    zcenters(i) = TransformStruct(i).Oc(3, 4);   
    
end

% Plot
if strcmp(plotoption, 'on') == 1
    plot3(xcenters(:, 1), ycenters(:, 1), zcenters(:, 1), 'Color', 'k', ...
        'Linewidth', 4)
end

% We can now use the list of TransformStruct.Oc to provide us information
% on the new centroids of all of these spheres (joints)

% New visualization
if strcmp(plotoption, 'on') == 1
    hold off
    figure()
    set(gcf, 'color', 'w')
    hold on
end

for i = 1:N+1  
    
    [TransformStruct(i).demonew] = SphericalSampling(TransformStruct(i).Oc(:, 4), ...
    TransformStruct(i).rs, colorvector(i, :), plotoption); 
    
    if strcmp(plotoption, 'on') == 1
        % Plot vectors which demonstrate the new Oc 
        quiver3(TransformStruct(i).Oc(1,4), TransformStruct(i).Oc(2,4), ...
            TransformStruct(i).Oc(3,4), TransformStruct(i).Oc(1,1), ...
            TransformStruct(i).Oc(2,1), TransformStruct(i).Oc(3,1), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'red');

        quiver3(TransformStruct(i).Oc(1,4), TransformStruct(i).Oc(2,4), ...
            TransformStruct(i).Oc(3,4), TransformStruct(i).Oc(1,2), ...
            TransformStruct(i).Oc(2,2), TransformStruct(i).Oc(3,2), ...
            'AutoScaleFactor', 0.05, 'Linewidth', 1.1, 'Color', 'green');

        quiver3(TransformStruct(i).Oc(1,4), TransformStruct(i).Oc(2,4), ...
            TransformStruct(i).Oc(3,4), TransformStruct(i).Oc(1,3), ...
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