% Kinegami Test V1 (Spherical Wrist)
% Last Edited 6/29/2021 by Lucien Peach

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

[infostruct, TransformStruct] = Kinegami(D, r, n, JointStruct, mirror, ...
    triple, theta_mod, fingertip);


