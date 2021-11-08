% Kinegami Test V1 (Cylindrical)
% Last Edited 7/20/2021 by Lucien Peach

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

% Specify whether DXF generation and save file should occur ('on'/'off')
DXF = 'on';

% Specify whether elbow splitting should occur past pi/2 ('on'/'off')
split = 'on';

% Input the kinematic chain robot specifications
% Specify DH Parameters, 
D = [0, 0.0001,    0.1,      0; ...
     0,      0,   0.08,      0; ...
     0,      0,    0.1,      0];
% Number of joints
n = size(D,1);

% Specify joint informations:
% Types of joints: 'R': Revolute joint, 'P': Prismatic joint, ...
% 'F': Fingertip, 'V': Vertex (not a joint)
TYPE = ['R', 'P', 'F']; 

% Maximum joint range
Qm = [pi, 0, pi];

% Initial joint configuration (last column of the DH table)
Q0 = [0, 0.04, 0];

% Specify the angle modification utilized (recommended: zeros(n))
theta_mod = [0, 0, 0];

% Layer of recursive sink gadget for revolute joint 
Nz = [1, 1, 1];

% Specify the orientation of the fingertip: 'x', 'y', or 'z'
fingertip = 'z';
 
% Specify number of sides for the polygon base of the prism tube
nsides = 4;

% Specify radius [m]
r = 0.02;

% Initialize JointStruct
JointStruct(n) = struct();
for i = 1:n
    JointStruct(i).qm = Qm(i);
    JointStruct(i).q0 = Q0(i);
    JointStruct(i).type = TYPE(i);
    JointStruct(i).nz = Nz(i);
end
N = size(JointStruct, 2) - 1;
TransformStruct(N+1) = struct();

% If the selfassign tag is applied, provide Oc for each joint
if strcmp(selfassign, 'true') == 1
    
    % Specify the orientation of the fingertip: 'x', 'y', or 'z'
    fingertip = 'x';
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
    % Otherwise, do nothing new
end

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, nsides, JointStruct, ...
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split);
