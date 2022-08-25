% KINEGAMI_HOPPER - Generates a hopper prototype with Kinegami.

% Authors: 
% Yiyu Wu <yiyuwu@sas.upenn.edu >
% Wei-Hsi Chen <weicc@seas.upenn.edu>
% Last Edited 1/1/2022
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


clear
close all
clc

% Add all the folders and subfolders to the search path
addpath(genpath(fileparts(mfilename('fullpath'))));

%% User Options - Change Prior to Running (if necessary)

% Determines which joint placement method will be implemented:
% For the general method, use JointPlacementA.m ('placementA'); 
% To guarantee no self-intersection, use JointPlacementB.m ('placementB');
% To locate the joints manually, use SelfAssign.m ('selfassign')
jointselect = 'selfassign';

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

% Specify whether or not pre-segmentation (for printing) is enabled
% ('on'/'off')
segmentation = 'off';

% Specify whether intermediary plots should be run ('on'/'off'). 'off' is
% recommended for faster computational time, 'on' is recommended for more
% in-depth analysis.
plotoption = 'off';

% Input the kinematic chain robot specifications
% Specify DH Parameters, 
a3 = 0.08;
d4 = 0.08;
d = 0.1;
D = [0,      0,      0,  -pi/2; ...
     0,  -pi/2,      0,      0; ...
%      a3,     0,      0,      0; ...
%      0,  -pi/2,     d4,      0; ...
%      0,   pi/2,      0,   pi/2; ...
%      0,   pi/2,      0,   pi/2; ...
     0,      0,   0.16,      0];
% Number of joints
n = size(D,1);

% Specify joint informations:
% Types of joints: 'R': Revolute joint, 'P': Prismatic joint, ...
% 'F': Fingertip, 'V': Vertex (not a joint)
TYPE = ['V', 'R', 'F']; 

% Maximum joint range
Qm = [pi, pi, pi];

% Initial joint configuration (last column of the DH table)
Q0 = D(:,end)';

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
if strcmp(jointselect, 'selfassign') == 1
    
    % Specify the orientation of the fingertip: 'x', 'y', or 'z'
    fingertip = 'x';
    
%     TransformStruct(1).Oc = [0, 1, 0, 0; ...
%                              0, 0, 1, 0; ...
%                              1, 0, 0, 0];
    
    TransformStruct(1).Oc = [1, 0, 0,  3/4*d; ...
        0, 1, 0, 0; ...
        0, 0, -1, -1/3*d];
    
%     TransformStruct(2).Oc = [-1/sqrt(2), 1/sqrt(2), 0, d; ...
%                               0, 0, -1, 0; ...
%                               1/sqrt(2), 1/sqrt(2), 0, 0];
    
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
    elbow_tuck, triple, theta_mod, fingertip, TransformStruct, ...
    DXF, split, segmentation, plotoption, jointselect);