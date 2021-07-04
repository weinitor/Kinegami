% Kinegami Test V1 (Planar)
% Last Edited 6/29/2021 by Lucien Peach

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

[infostruct, TransformStruct] = Kinegami(D, r, n, JointStruct, mirror, ...
    triple, theta_mod, fingertip);
