% JOINTTEST - Used for running and intermediary visualization of
% JointPlacement.m.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last Edited 1/17/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


clear
close all
clc

% Add all the folders and subfolders to the search path
addpath(genpath(fileparts(mfilename('fullpath'))));

% User Options - Change Prior to Running (if necessary)

% Determines whether the user wishes to use DH parameters ('false') or
% assign the Joint Parameters themselves ('true')
selfassign = 'false';

% Determines whether the user wishes to have elbow joints mirrored ('on')
% or appear normally ('off')
elbow_tuck = 'on';

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
split = 'on';

% Specify whether or not pre-segmentation (for printing) is enabled
% ('on'/'off')
segmentation = 'on';

% Specify whether intermediary plots should be run ('on'/'off'). 'off' is
% recommended for faster computational time, 'on' is recommended for more
% in-depth analysis.
plotoption = 'on';

% Specify DH Parameters, if needed
D = [0.15, 0, 0, 0; ...
    0.15, 0, 0, 0; ...
    0.15, 0, 0, 0; ...
    0.15, 0, 0, 0];

% Specify number of sides (polygon)
nsides = 4;

% Specify radius [m]
r = 0.02;

% Specify number of joints
n = 4;

JointStruct(4) = struct();

for i = 1:4
    JointStruct(i).qm = 4/3*pi;
    JointStruct(i).q0 = 0;
    JointStruct(i).type = 'R';
    JointStruct(i).nz = 1;
end

JointStruct(4).qm = pi;
JointStruct(4).type = 'F';

N = size(JointStruct, 2) - 1;

% Otherwise, do nothing besides initialization
TransformStruct(N+1) = struct();

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, n, JointStruct, ...
    elbow_tuck, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split, segmentation, plotoption);

