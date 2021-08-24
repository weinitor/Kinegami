% Kinegami Test V1 (Cylindrical)
% Last Edited 7/20/2021 by Lucien Peach

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
fingertip = 'z';

% Specify whether DXF generation and save file should occur ('on'/'off')
DXF = 'on';

% Specify whether elbow splitting should occur past pi/2 ('on'/'off')
split = 'on';

% Specify DH Parameters, if needed
D = [0, 0.0001, 0.1, 0; ...
    0, 0, 0.08, 0; ...
    0, 0, 0.1, 0];

% Specify radius [m]
r = 0.02;

% Specify number of joints
n = 4;

JointStruct(1).q0 = 0;
JointStruct(1).qm = pi;
JointStruct(1).type = 'R';

JointStruct(2).q0 = 0.04;
JointStruct(2).qm = 0;
JointStruct(2).type = 'P';

JointStruct(3).q0 = 0;
JointStruct(3).qm = pi;
JointStruct(3).type = 'F';

N = size(JointStruct, 2) - 1;

% If the selfassign tag is applied, provide Oc for each joint
if strcmp(selfassign, 'true') == 1
    
    TransformStruct(N+1) = struct();
    
    TransformStruct(1).Oc = [1, 0, 0, 0; ...
        0, 1, -0.01, 0; ...
        0, 0.01, 1, 0.0176];
    
    TransformStruct(2).Oc = [0, 1, 0, 0; ...
        -0.01, 0, 1, 0; ...
        1, 0, 0.01, 0.1459];
    
    TransformStruct(3).Oc = [0, 1, 0, 0; ...
        -0.01, 0, 1, 0; ...
        1, 0, 0.01, 0.28];    
    
else
    
    % Otherwise, do nothing besides initialization
    TransformStruct(N+1) = struct();

end

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, n, JointStruct, ...
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split);

