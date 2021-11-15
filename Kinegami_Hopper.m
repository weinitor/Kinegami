% Kinegami Test V1 (6 DOF)
% Last Edited 11/5/2021 by Shelly Wu and Wei-Hsi Chen

clear
close all
clc

%% USER OPTIONS - Change Prior to Running (if necessary)

% Determines whether the user wishes to use DH parameters ('false') or
% assign the Joint Parameters themselves ('true')
selfassign = 'true';

% Determines whether the user wishes to have elbow joints mirrored ('on')
% or appear normally ('off')
mirror = 'on';

% Specify whether elbow splitting should occur past pi/2 ('on'/'off')
split = 'on';

% Determines whether the user wishes to print 3 iterations of the print
% pattern ('triple' - recommended) or 2 ('double')
triple = 'triple';

% Specify whether DXF generation and save file should occur ('on'/'off')
DXF = 'on';

% Specify whether or not pre-segmentation (for printing) is enabled
% ('on'/'off')
segmentation = 'off';

%% KINEMATIC CHAIN SPECIFICATION

% Specify number of sides and the circumradius [m] for the polygon base 
% of the prism tube
nsides = 4;
r = 0.02;

% Input the kinematic chain robot specifications
% Number of joints and initialize DH Parameter table D
% n equals to the number of joints plus one fingertip
n = 3;
D = zeros(n,4);

% Specify DH Parameters nX4, in the order of "Link length (a)", 
% "Link twist (alpha)", "Joint offset (d)", and "Joint angle (theta)".
d = 0.1;
D = [0,      0,      0,  -pi/2; ...
     0,  -pi/2,      0,      0; ...
     0,      0,   0.16,      0];
 
% Specify joint information as a row vector:
% Contains n elements for each vector
% Types of joints: 'R': Revolute joint, 'P': Prismatic joint, ...
% 'F': Fingertip, 'V': Spine Vertex (not a joint)
TYPE = ['V', 'R', 'F']; 

% Maximum joint range (row vec.)
Qm = [pi, pi, pi];

% Initial joint configuration (last column of the DH table)
Q0 = D(:,end);

% Specify the angle modification utilized (recommended: zeros(n)) (row vec.)
theta_mod = [0, 0, 0];

% Layer of recursive sink gadget for revolute joint (row vec.)
Nz = [1, 1, 1];

% Specify the orientation of the fingertip: 'x', 'y', or 'z'
fingertip = 'z';

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

%% SELF-ASSIGNED JOINT LOCATION

% If the selfassign tag is applied, provide Oc for each joint
if strcmp(selfassign, 'true') == 1
    
    % Specify the orientation of the fingertip: 'x', 'y', or 'z'
    fingertip = 'x';
    
    TransformStruct(1).Oc = [1, 0, 0,  3/4*d; ...
        0, 1, 0, 0; ...
        0, 0, -1, -1/3*d];
    
    TransformStruct(2).Oc = [-1/sqrt(2), 1/sqrt(2), 0, 3/4*d; ...
        0, 0, -1, 0; ...
        1/sqrt(2), 1/sqrt(2), 0, 1/4*d];  
    
    TransformStruct(3).Oc = [-1/sqrt(2), 1/sqrt(2), 0, 1/3*d; ...
        0, 0, -1, 0; ...
        1/sqrt(2), 1/sqrt(2), 0, 2/3*d];
else  
    % Otherwise, do nothing new
end

%% RUN KINEGAMI

% Run Kinegami code
[infostruct, TransformStruct, DataNet] = Kinegami(D, r, nsides, JointStruct, ...
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split, segmentation);