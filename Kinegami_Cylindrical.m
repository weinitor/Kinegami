% Kinegami Test V1 (Cylindrical)
% Last Edited 7/20/2021 by Lucien Peach

clear
close all
clc

D = [0, 0.0001, 0.1, 0; ...
    0, 0, 0.08, 0; ...
    0, 0, 0.1, 0];

r = 0.02;
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

mirror = 'on';
triple = 'triple';
theta_mod = [0, 0, 0, 0];
fingertip = 'z';
selfassign = 'true';

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

DXF = 'on';
split = 'on';

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, n, JointStruct, ...
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split);

