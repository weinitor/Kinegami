% KINEGAMI_KNOT_TEST - This example is to show the possibility of
% self-intersection. The comments for the joint specifications are deleted
% for simplicity. 

% Authors: 
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Last Edited 1/26/2022
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.txt for detail.


clear
close all
clc

%% USER OPTIONS - Change Prior to Running (if necessary)

% Determines which joint placement method will be implemented:
% For the general method, use JointPlacementA.m ('placementA'); 
% To guarantee no self-intersection, use JointPlacementB.m ('placementB');
% To locate the joints manually, use SelfAssign.m ('selfassign')
jointselect = 'placementB';

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
n = 3;
D = [12*0.02,      0,      0,      0; ...
     10*0.02,  -pi/2,      0,     pi; ...
     0,           0,  6*0.02,     0];
fingertip = 'z';
TYPE = ['R', 'R', 'F']; 
Qm = [pi, pi, pi]; 
Q0 = [0, 0, 0]; 
theta_mod = [0, 0, 0];
Nz = [1, 1, 1]; 

% % Universal joint
% n = 4;
% D = [0,           0,    8*0.02,          0; ...
%      0,        pi/2,         0,      -pi/2; ...
%      0,       -pi/2,         0,       0*pi; ...
%      -6*0.02,  pi/2,         0,        pi];
% fingertip = 'x';
% TYPE = ['R', 'R', 'R', 'F']; 
% Qm = [pi, pi, pi, pi]; 
% Q0 = [0, 0, 0, 0]; 
% theta_mod = [pi/2, 0, 0, 0]; 
% Nz = [1, 1, 1, 1]; 


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
if strcmp(jointselect, 'selfassign') == 1
    
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
[infostruct, TransformStruct, DataNet, JointStruct] = Kinegami(D, r, nsides, JointStruct, ...
    elbow_tuck, triple, theta_mod, fingertip, TransformStruct, ...
    DXF, split, segmentation, plotoption, jointselect);
