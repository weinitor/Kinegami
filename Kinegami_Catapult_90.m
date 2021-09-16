% Catapult V1.2 (theta_m = 90deg)
% Last Edited 9/8/2021 by Lucien Peach

clear
close all
clc

% User Options - Change Prior to Running (if necessary)

% Determines whether the user wishes to use DH parameters ('false') or
% assign the Joint Parameters themselves ('true')
selfassign = 'true';

% Determines whether the user wishes to have elbow joints mirrored ('on')
% or appear normally ('off')
mirror = 'on';

% Determines whether the user wishes to print 3 iterations of the print
% pattern ('triple' - recommended) or 2 ('double')
triple = 'triple';

% Specify the angle modification utilized ([0, 0, 0, 0] recommended)
theta_mod = [0, 0, 0, 0];

% Specify the orientation of the fingertip: 'x', 'y', or 'z'
fingertip = 'x';

% Specify whether DXF generation and save file should occur ('on'/'off')
DXF = 'on';

% Specify whether elbow splitting should occur past pi/2 ('on'/'off')
split = 'off';

% Specify number of sides (polygon)
nsides = 4;

% Specify radius [m]
r = 0.05;

% Specify DH Parameters, if needed
D = [0, pi/2, 0, pi/2; ...
    0, pi/2, 0, 0; ...
    0.2, 0, 0, 0];

% Specify number of joints
n = 6;

if strcmp(selfassign, 'false') == 1

    n = 3;

    JointStruct(n) = struct();

    for i = 1:n
        JointStruct(i).qm = pi;
        JointStruct(i).q0 = 0;
        JointStruct(i).type = 'R';
    end

    JointStruct(1).q0 = pi/2;
    JointStruct(1).type = 'F';

    N = size(JointStruct, 2) - 1;
    
end

% If the selfassign tag is applied, provide Oc for each joint
% Make sure that Oc(:,4) are all not equal to 0 (x.xxx * 10^-25, etc., is
% acceptable)
if strcmp(selfassign, 'true') == 1
    
    d = 0.1;
    
    JointStruct(n) = struct();

    for i = 1:n
        JointStruct(i).qm = pi/2;
        JointStruct(i).q0 = 0;
        JointStruct(i).type = 'V';
    end
    
    JointStruct(5).type = 'R';
    JointStruct(5).nz = 2;
    
    JointStruct(6).type = 'F';
    

    N = size(JointStruct, 2) - 1;
    
    TransformStruct(N+1) = struct();
    
    TransformStruct(1).Oc = [1, 0, 0, 0; ...
        0, 0, -1, d; ...
        0, 1, 0, 0];
    
    TransformStruct(2).Oc = [0, 0, 1, 2*d; ...
        1, 0, 0, 0; ...
        0, 1, 0, 0];
    
    TransformStruct(3).Oc = [-1, 0, 0, 0; ...
        0, 0, 1, -d; ...
        0, 1, 0, 0];  
    
    TransformStruct(4).Oc = [1, 0, 0, 0; ...
        0, 0, -1, 0; ...
        0, 1, 0, 0];
    
    TransformStruct(5).Oc = [0, -1, 0, 0.5*d; ...
        0, 0, -1, 0; ...
        1, 0, 0, d];

    TransformStruct(6).Oc = [0, -1, 0, 0.51*d; ...
        0, 0, -1, 0; ...
        1, 0, 0, 3*d];
    
    % Note in lines 112-114 an important specification. Due to singularity
    % issues, TransformStruct(6).Oc(1,4) = 0.51*d instead of 0.5*d.
    % Referencing the 2D Crease Pattern, we can see that this minor change
    % creates an elbow joint outside of our tolerance range, and thus only
    % a tube is generated. Therefore, the functionality of the final
    % origami folded model is not impacted by this necessary change. 
    
    % Here solely to not mess up Kinegami running
    D = [0, pi/2, 0, pi/2; ...
        0, pi/2, 0, 0; ...
        0.2, 0, 0, 0];
    
else
    
    % Otherwise, do nothing besides initialization
    TransformStruct(N+1) = struct();

end

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, nsides, JointStruct, ...
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split);
