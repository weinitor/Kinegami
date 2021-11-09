% Kinegami Test V1 (Planar)
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

% If the selfassign tag is applied, provide Oc for each joint
% Make sure that Oc(:,4) are all not equal to 0 (x.xxx * 10^-25, etc., is
% acceptable)
if strcmp(selfassign, 'true') == 1
    
    TransformStruct(N+1) = struct();
    
    TransformStruct(1).Oc = [1, 0, 0, 0.15; ...
        0, 1, 0, 0; ...
        0, 0, 1, 7.0842*10^-9];
    
    TransformStruct(2).Oc = [-1, 0, 0, 0.3; ...
        0, 1, 0, 0; ...
        0, 0, -1, -4.5539*10^-10];
    
    TransformStruct(3).Oc = [1, 0, 0, 0.45; ...
        0, 1, 0, 0; ...
        0, 0, 1, -6.2523*10^-9]; 
    
    TransformStruct(4).Oc = [1, 0, 0, 0.6; ...
        0, 1, 0, 0; ...
        0, 0, 1, 0];
    
else
    
    % Otherwise, do nothing besides initialization
    TransformStruct(N+1) = struct();

end

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, n, JointStruct, ...
    elbow_tuck, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split, segmentation);
