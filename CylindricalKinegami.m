% Kinegami Test V1 (Cylindrical)
% Last Edited 7/1/2021 by Lucien Peach

clear
close all
clc

D = [0, 0.01, 0.1, 0; ...
    0, 0, 0.1, 0; ...
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
JointStruct(3).qm = 0;
JointStruct(3).type = 0;


mirror = 'on';
triple = 'triple';
theta_mod = [0, 0, 0, 0];
fingertip = 'z';

[infostruct, TransformStruct] = Kinegami(D, r, n, JointStruct, mirror, ...
    triple, theta_mod, fingertip);

