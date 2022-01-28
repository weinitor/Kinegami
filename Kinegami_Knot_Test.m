% Kinegami Test V1 (Intersection Testing)
% Last Edited 1/26/2022 by Wei-Hsi
% This example is to show the possibility of self-intersection

clear
close all
clc

%% USER OPTIONS - Change Prior to Running (if necessary)

% Determines whether the user wishes to use DH parameters ('false') or
% assign the Joint Parameters themselves ('true')
selfassign = 'false';

% Determines whether the user wishes their elbow fittings to have visible
% tucks ('on' - recommended) or appear with only the lower outlines ('off')
elbow_tuck = 'on';

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

% Specify whether intermediary plots should be run ('on'/'off'). 'off' is
% recommended for faster computational time, 'on' is recommended for more
% in-depth analysis.
plotoption = 'off';

%% KINEMATIC CHAIN SPECIFICATION (Universal)

% Specify number of sides and the circumradius [m] for the polygon base 
% of the prism tube
nsides = 4;
r = 0.02;

% Input the kinematic chain robot specifications
% Number of joints and initialize DH Parameter table D
% n equals to the number of joints plus one fingertip

% Specify DH Parameters nX4, in the order of "Link length (a)", 
% "Link twist (alpha)", "Joint offset (d)", and "Joint angle (theta)".

% 2 R arm
% n = 3;
% D = [12*0.02,      0,      0,      0; ...
%      10*0.02,  -pi/2,      0,     pi; ...
%      0,           0,  6*0.02,     0];
% fingertip = 'z';
% TYPE = ['R', 'R', 'F']; 
% Qm = [pi, pi, pi]; 
% Q0 = [0, 0, 0]; 
% theta_mod = [0, 0, 0];
% Nz = [1, 1, 1]; 

% % Universal joint
n = 4;
D = [0,           0,    8*0.02,          0; ...
     0,        pi/2,         0,      -pi/2; ...
     0,       -pi/2,         0,       0*pi; ...
     -6*0.02,  pi/2,         0,        pi];
fingertip = 'x';

% Specify joint information as a row vector:
% Contains n elements for each vector
% Types of joints: 'R': Revolute joint, 'P': Prismatic joint, ...
% 'F': Fingertip, 'V': Spine Vertex (not a joint)
TYPE = ['V', 'R', 'R', 'F']; 

% Maximum joint range (row vec.)
Qm = [pi, pi, pi, pi]; 

% Initial joint configuration (row vec.)
Q0 = [0, 0, 0, 0]; 

% Specify the angle modification utilized (recommended: zeros(n)) (row vec.)
theta_mod = [pi/2, 0, 0, 0]; 

% Layer of recursive sink gadget for revolute joint (row vec.)
Nz = [1, 1, 1, 1]; 

% Specify the orientation of the fingertip: 'x', 'y', or 'z'
% fingertip = 'z';

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
    
    n = 4;
    TYPE = ['V', 'R', 'R', 'F']; 
    Qm = [pi, pi, pi, pi]; 
    Q0 = [0, pi/2, 0, 0]; 
    theta_mod = [0, 0, 0, 0]; 
    Nz = [1, 1, 1, 1]; 
    JointStruct = struct() % reset struct
    JointStruct(n) = struct();
    for i = 1:n
        JointStruct(i).qm = Qm(i);
        JointStruct(i).q0 = Q0(i);
        JointStruct(i).type = TYPE(i);
        JointStruct(i).nz = Nz(i);
    end
    N = size(JointStruct, 2) - 1;
    TransformStruct(N+1) = struct();
    
    % Specify the orientation of the fingertip: 'x', 'y', or 'z'
    fingertip = 'x';
    
    % Specify the frame (3X4) of each individul joint
    TransformStruct(1).Oc = [0, 0, 1, 0; ...
        0, -1, 0, 0; ...
        1, 0, 0, -4*r];
    
    TransformStruct(2).Oc = [0, 0, 1, -2.5*r; ...
        1, 0, 0, 0; ...
        0, 1, 0, 0];
    
    TransformStruct(3).Oc = [1, 0, 0, 0; ...
        0, 0, -1, 2.5*r; ...
        0, 1, 0, 0];  
    
    TransformStruct(4).Oc = [0, 0, 1, 0; ...
        0, -1, 0, 0; ...
        1, 0, 0, 4*r];            
    % etc...
    
else  
    % Otherwise, do nothing new
end

%% RUN KINEGAMI

% Run Kinegami code
[infostruct, TransformStruct, DataNet] = Kinegami(D, r, nsides, JointStruct, ...
    elbow_tuck, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split, segmentation, plotoption);
