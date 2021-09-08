% Kinegami Test V1 (Spherical Wrist)
% Last Edited 7/21/2021 by Lucien Peach

clear
close all
clc

% User Options - Change Prior to Running (if necessary)

% Determines whether the user wishes to use DH parameters ('false') or
% assign the Joint Parameters themselves ('true')
selfassign = 'false';

% Determines whether the user wishes to have elbow joints mirrored ('on')
% or appear normally ('off')
mirror = 'on';

% Determines whether the user wishes to print 3 iterations of the print
% pattern ('triple' - recommended) or 2 ('double')
triple = 'triple';

% Specify the angle modification utilized ([0, 0, 0, 0] recommended)
theta_mod = [0, 0, 0, 0];

% Specify the orientation of the fingertip: 'x', 'y', or 'z'
fingertip = 'z';

% Specify whether DXF generation and save file should occur ('on'/'off')
DXF = 'on';

% Specify whether elbow splitting should occur past pi/2 ('on'/'off')
split = 'on';

% Specify DH Parameters, if needed
D = [0, 0, 0.1, 0; ...
    0, pi/2, 0, pi/2; ...
    0, pi/2, 0, pi/2; ...
    0, pi/2, 0.1, 0];

% Specify number of sides (polygon)
nsides = 4;

% Specify radius [m]
r = 0.02;

% Specify number of sides of polygon
n = 4;

JointStruct(4) = struct();

for i = 1:4
    JointStruct(i).qm = pi;
    JointStruct(i).q0 = 0;
    JointStruct(i).type = 'R';
    JointStruct(i).nz = 1;
end

JointStruct(2).q0 = pi/2;
JointStruct(3).q0 = pi/2;
JointStruct(4).type = 'F';

N = size(JointStruct, 2) - 1;

% If the selfassign tag is applied, provide Oc for each joint
if strcmp(selfassign, 'true') == 1
    
    TransformStruct(N+1) = struct();
    
    TransformStruct(1).Oc = [1, 0, 0, 0; ...
        0, 1, 0, 0; ...
        0, 0, 1, -0.159];
    
    TransformStruct(2).Oc = [0, -1, 0, 0; ...
        0, 0, -1, 0.2377; ...
        1, 0, 0, 0.1];
    
    TransformStruct(3).Oc = [0, 0, 1, -0.139; ...
        -1, 0, 0, 0; ...
        0, -1, 0, 0.1];    
    
    TransformStruct(4).Oc = [0, 0, 1, 0; ...
        0, -1, 0, 0; ...
        1, 0, 0, 0.2];
    
else
    
    % Otherwise, do nothing besides initialization
    TransformStruct(N+1) = struct();

end

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, n, JointStruct, ...
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split);


