% Kinegami Test V1 (6 DOF)
% Last Edited 9/9/2021 by Wei-Hsi Chen

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

% Specify whether DXF generation and save file should occur ('on'/'off')
DXF = 'on';

% Specify whether elbow splitting should occur past pi/2 ('on'/'off')
split = 'on';

% Input the kinematic chain robot specifications
% Specify DH Parameters, 
a3 = 0.08;
d4 = 0.08;
D = [0,      0,      0,  -pi/2; ...
     0,  -pi/2,      0,      0; ...
     a3,     0,      0,      0; ...
     0,  -pi/2,     d4,      0; ...
     0,   pi/2,      0,   pi/2; ...
     0,   pi/2,      0,   pi/2; ...
     0,      0,   0.16,      0];
% Number of joints
n = size(D,1);

% Specify joint informations:
% Types of joints: 'R': Revolute joint, 'P': Prismatic joint, ...
% 'F': Fingertip, 'V': Vertex (not a joint)
TYPE = ['R', 'R', 'R', 'R', 'R', 'R', 'F']; 

% Maximum joint range
Qm = [pi, pi, pi, pi, pi, pi, pi];

% Initial joint configuration (last column of the DH table)
Q0 = D(:,end);

% Specify the angle modification utilized (recommended: zeros(n))
theta_mod = [0, 0, 0, 0, 0, 0, 0];

% Layer of recursive sink gadget for revolute joint 
Nz = [1, 1, 1, 1, 1, 1, 1];

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
    
%     TransformStruct(1).Oc = [0, 1, 0, 0; ...
%                              0, 0, 1, 0; ...
%                              1, 0, 0, 0];
    
    TransformStruct(1).Oc = [-1, 0, 0, 0; ...
                             0, -1, 0, 0; ...
                             0, 0, -1, 6*r];
    
    TransformStruct(2).Oc = [1, 0, 0, 0; ...
                             0, 0, 1, 0; ...
                             0, -1, 0, 12*r];  
    
    TransformStruct(3).Oc = [0, -1, 0, a3; ...
                             0, 0, 1, 0; ...
                             -1, 0, 0, 12*r];
    
    TransformStruct(4).Oc = [0, -1, 0, a3; ...
                             -1, 0, 0, 0; ...
                             0, 0, -1, 12*r-d4];
    
    TransformStruct(5).Oc = [0, -1, 0, a3; ...
                             0, 0, 1, -4*r; ...
                             -1, 0, 0, 8*r-d4];  
    
    TransformStruct(6).Oc = [0, 0, 1, a3+4*r; ...
                             1, 0, 0, 0; ...
                             0, 1, 0, 8*r-d4];
                         
    TransformStruct(7).Oc = [1, 0, 0, a3+8*r; ...
                             0, 1, 0, 0; ...
                             0, 0, 1, 8*r-d4];
else  
    % Otherwise, do nothing new
end

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, nsides, JointStruct, ...
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split);
