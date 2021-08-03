% Kinegami Test V1 (Universal)
% Last Edited 7/21/2021 by Lucien Peach

clear
close all
clc

D = [0, pi/2, 0, pi/2; ...
    0, pi/2, 0, 0; ...
    0.2, 0, 0, 0];

r = 0.02;
n = 4;

JointStruct(3) = struct();

for i = 1:3
    JointStruct(i).qm = pi;
    JointStruct(i).q0 = 0;
    JointStruct(i).type = 'R';
end

JointStruct(1).q0 = pi/2;
JointStruct(3).type = 0;

mirror = 'on';
triple = 'triple';
theta_mod = [0, 0, 0, 0];
fingertip = 'x';
selfassign = 'false';

N = size(JointStruct, 2) - 1;

% If the selfassign tag is applied, provide Oc for each joint
% Make sure that Oc(:,4) are all not equal to 0 (x.xxx * 10^-25, etc., is
% acceptable)
if strcmp(selfassign, 'true') == 1
    
    TransformStruct(N+1) = struct();
    
    TransformStruct(1).Oc = [0, -1, 0, 0; ...
        0, 0, -1, 0.1843; ...
        1, 0, 0, -1.1297*10^-17];
    
    TransformStruct(2).Oc = [0, 0, 1, 6.9248*10^-9; ...
        0, -1, 0, -4.2402*10^-25; ...
        1, 0, 0, -4.2402*10^-25];
    
    TransformStruct(3).Oc = [0, 0, 1, 1.12246*10^-17; ...
        0, -1, 0, 1.12246*10^-17; ...
        1, 0, 0, 0.2];    
    
else
    
    % Otherwise, do nothing besides initialization
    TransformStruct(N+1) = struct();

end

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, n, JointStruct, ...
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct);

