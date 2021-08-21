% Kinegami Test V1 (Universal)
% Last Edited 7/21/2021 by Lucien Peach

clear
close all
clc

selfassign = 'true';

if strcmp(selfassign, 'false') == 1

    D = [0, pi/2, 0, pi/2; ...
        0, pi/2, 0, 0; ...
        0.2, 0, 0, 0];

    r = 0.02;
    n = 3;

    JointStruct(n) = struct();

    for i = 1:n
        JointStruct(i).qm = pi;
        JointStruct(i).q0 = 0;
        JointStruct(i).type = 'R';
    end

    JointStruct(1).q0 = pi/2;
    JointStruct(3).type = 'F';

    mirror = 'on';
    triple = 'triple';
    theta_mod = [0, 0, 0, 0];
    fingertip = 'x';

    N = size(JointStruct, 2) - 1;
    
end

% If the selfassign tag is applied, provide Oc for each joint
% Make sure that Oc(:,4) are all not equal to 0 (x.xxx * 10^-25, etc., is
% acceptable)
if strcmp(selfassign, 'true') == 1
    
    r = 0.02;
    n = 4;

    JointStruct(n) = struct();

    for i = 1:n
        JointStruct(i).qm = pi;
        JointStruct(i).q0 = 0;
        JointStruct(i).type = 'R';
    end

    JointStruct(2).q0 = pi/2;
    JointStruct(1).type = 'V';
    JointStruct(4).type = 'F';

    mirror = 'on';
    triple = 'triple';
    theta_mod = [0, 0, 0, 0];
    fingertip = 'x';

    N = size(JointStruct, 2) - 1;
    
    TransformStruct(N+1) = struct();
    
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
    
    % Here solely to not mess up Kinegami running
    D = [0, pi/2, 0, pi/2; ...
        0, pi/2, 0, 0; ...
        0.2, 0, 0, 0];
    
else
    
    % Otherwise, do nothing besides initialization
    TransformStruct(N+1) = struct();

end

DXF = 'on';
split = 'on';

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, n, JointStruct, ...
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split);



