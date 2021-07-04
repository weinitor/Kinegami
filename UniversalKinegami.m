% Kinegami Test V1 (Universal)
% Last Edited 6/29/2021 by Lucien Peach

clear
close all
clc

D = [0, pi/2, 0, -pi/2; ...
    0, pi/2, 0, pi/2; ...
    0, pi/2, 0.2, 0];

r = 0.02;
n = 4;

JointStruct(3) = struct();

for i = 1:3
    JointStruct(i).qm = pi;
    JointStruct(i).q0 = 0;
    JointStruct(i).type = 'R';
end

JointStruct(1).q0 = -pi/2;
JointStruct(2).q0 = pi/2;
JointStruct(3).type = 0;

mirror = 'on';
triple = 'triple';
theta_mod = [0, 0, 0, 0];
fingertip = 'z';

[infostruct, TransformStruct] = Kinegami(D, r, n, JointStruct, mirror, ...
    triple, theta_mod, fingertip);

