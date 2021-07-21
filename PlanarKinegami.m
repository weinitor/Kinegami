% Kinegami Test V1 (Planar)
% Last Edited 7/21/2021 by Lucien Peach

clear
close all
clc

D = [0.15, 0, 0, 0; ...
    0.15, 0, 0, 0; ...
    0.15, 0, 0, 0; ...
    0.15, 0, 0, 0];

r = 0.02;
n = 4;

JointStruct(4) = struct();

for i = 1:4
    JointStruct(i).qm = 4/3*pi;
    JointStruct(i).q0 = 0;
    JointStruct(i).type = 'R';
end

JointStruct(4).type = 0;

mirror = 'on';
triple = 'triple';
theta_mod = [0, 0, 0, 0];
fingertip = 'x';
selfassign = 'true';

N = size(JointStruct, 2) - 1;

% If the selfassign tag is applied, provide Oc for each joint
% Make sure that Oc(:,4) are all not equal to 0 (x.xxx * 10^-25, etc., is
% acceptable)
if strcmp(selfassign, 'true') == 1
    
    TransformStruct(N+1) = struct();
    
    TransformStruct(1).Oc = [1, 0, 0, 0.15; ...
        0, 1, 0, 0; ...
        0, 0, 1, 7.0842*10^-9];
    
    TransformStruct(2).Oc = [1, 0, 0, 0.3; ...
        0, 1, 0, 0; ...
        0, 0, 1, -4.5539*10^-10];
    
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
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct);
