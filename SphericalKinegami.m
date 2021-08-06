% Kinegami Test V1 (Spherical Wrist)
% Last Edited 7/21/2021 by Lucien Peach

clear
close all
clc

D = [0, 0, 0.1, 0; ...
    0, pi/2, 0, pi/2; ...
    0, pi/2, 0, pi/2; ...
    0, pi/2, 0.1, 0];

r = 0.03;
n = 4;

JointStruct(4) = struct();

for i = 1:4
    JointStruct(i).qm = pi;
    JointStruct(i).q0 = 0;
    JointStruct(i).type = 'R';
end

JointStruct(2).q0 = pi/2;
JointStruct(3).q0 = pi/2;
JointStruct(4).type = 0;

mirror = 'on';
triple = 'triple';
theta_mod = [0, 0, 0, 0];
fingertip = 'z';
selfassign = 'false';

N = size(JointStruct, 2) - 1;

% If the selfassign tag is applied, provide Oc for each joint
if strcmp(selfassign, 'true') == 1
    
    TransformStruct(N+1) = struct();
    
    TransformStruct(1).Oc = [1, 0, 0, 0; ...
        0, 1, 0, 0; ...
        0, 0, 1, -0.159];
    
    TransformStruct(2).Oc = [0, -1, 0, 0; ...
        0, 0, -1, 0.2377; ...
        1, 0, 0, 0.1];
    
    TransformStruct(3).Oc = [0, 0, 1, -0.139; ...
        -1, 0, 0, 0; ...
        0, -1, 0, 0.1];    
    
    TransformStruct(4).Oc = [0, 0, 1, 0; ...
        0, -1, 0, 0; ...
        1, 0, 0, 0.2];
    
else
    
    % Otherwise, do nothing besides initialization
    TransformStruct(N+1) = struct();

end

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, n, JointStruct, ...
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct);


