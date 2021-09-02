% Kinegami Test V1 (6 DOF)
% Last Edited 7/21/2021 by Wei-Hsi Chen

clear
close all
clc

% User Options - Change Prior to Running (if necessary)

% Determines whether the user wishes to use DH parameters ('false') or
% assign the Joint Parameters themselves ('true')
selfassign = 'false';

% Determines whether the user wishes to have elbow joints mirrored ('on')
% or appear normally ('off')
mirror = 'on';

% Determines whether the user wishes to print 3 iterations of the print
% pattern ('triple' - recommended) or 2 ('double')
triple = 'triple';

% Specify the angle modification utilized ([0, 0, 0, 0] recommended)
theta_mod = [0, 0, 0, 0, 0, 0, 0];

% Specify the orientation of the fingertip: 'x', 'y', or 'z'
fingertip = 'x';

% Specify whether DXF generation and save file should occur ('on'/'off')
DXF = 'on';

% Specify whether elbow splitting should occur past pi/2 ('on'/'off')
split = 'on';

% Specify DH Parameters, if needed
a3 = 0.2;
d4 = 0.2;
D = [0,     0,  0, 0; ...
     0, -pi/2,  0, 0; ...
     a3,    0,  0, 0; ...
     0, -pi/2, d4, 0; ...
     0,  pi/2,  0, pi/2; ...
     0, pi/2,  0, pi/2; ...
     0,     0,  0, 0];
 
% Specify number of sides (polygon)
nsides = 4;

% Specify radius [m]
r = 0.02;

% Specify number of joints
n = 7;

if strcmp(selfassign, 'false') == 1

    JointStruct(n) = struct();

    for i = 1:n
        JointStruct(i).qm = pi;
        JointStruct(i).q0 = 0;
        JointStruct(i).type = 'R';
    end

    JointStruct(5).q0 = pi/2;
    JointStruct(6).q0 = pi/2;
    JointStruct(n).type = 'F';

    N = size(JointStruct, 2) - 1;
    
end

% If the selfassign tag is applied, provide Oc for each joint
% Make sure that Oc(:,4) are all not equal to 0 (x.xxx * 10^-25, etc., is
% acceptable)
if strcmp(selfassign, 'true') == 1
    
    JointStruct(n) = struct();

    for i = 1:n
        JointStruct(i).qm = pi;
        JointStruct(i).q0 = 0;
        JointStruct(i).type = 'R';
    end

    JointStruct(7).type = 'F';

    N = size(JointStruct, 2) - 1;
    
    TransformStruct(N+1) = struct();
    
%     TransformStruct(1).Oc = [0, 1, 0, 0; ...
%                              0, 0, 1, 0; ...
%                              1, 0, 0, 0];
    
    TransformStruct(1).Oc = [-1, 0, 0, 0; ...
                             0, -1, 0, 0; ...
                             0, 0, -1, 6*r];
    
    TransformStruct(2).Oc = [1, 0, 0, 0; ...
                             0, 0, 1, 0; ...
                             0, -1, 0, 12*r];  
    
    TransformStruct(3).Oc = [0, -1, 0, 4*r; ...
                             0, 0, 1, 0; ...
                             -1, 0, 0, 12*r];
    
    TransformStruct(4).Oc = [0, -1, 0, 4*r; ...
                             -1, 0, 0, 0; ...
                             0, 0, -1, 8*r];
    
    TransformStruct(5).Oc = [0, -1, 0, 4*r; ...
                             0, 0, 1, -4*r; ...
                             -1, 0, 0, 4*r];  
    
    TransformStruct(6).Oc = [0, 0, 1, 8*r; ...
                             1, 0, 0, 0; ...
                             0, 1, 0, 4*r];
                         
    TransformStruct(7).Oc = [1, 0, 0, 12*r; ...
                             0, 1, 0, 0; ...
                             0, 0, 1, 4*r];
                        
    
    % Here solely to not mess up Kinegami running
%     D = [0, pi/2, 0, pi/2; ...
%         0, pi/2, 0, 0; ...
%         0.2, 0, 0, 0];
D = [0,     0,  0, 0; ...
     0, -pi/2,  0, 0; ...
     a3,    0,  0, 0; ...
     0, -pi/2, d4, 0; ...
     0,  pi/2,  0, pi/2; ...
     0, pi/2,  0, pi/2; ...
     0,     0,  0, 0];
 
else
    
    % Otherwise, do nothing besides initialization
    TransformStruct(N+1) = struct();

end

[infostruct, TransformStruct, DataNet] = Kinegami(D, r, nsides, JointStruct, ...
    mirror, triple, theta_mod, fingertip, selfassign, TransformStruct, ...
    DXF, split);

