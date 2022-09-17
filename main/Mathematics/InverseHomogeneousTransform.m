function [DataStruct] = InverseHomogeneousTransform(index, D)
% INVERSEHOMOGENEOUSTRANSFORM - Finds the inverse transformation matrix 
% given the D-H specification.

% Inputs:
%   D           - the D-H parameter table of values: i x [a, alpha, d, 
%                 theta].
%   index       - ranges from 1 to i, used for identifying and properly 
%                 placing T into a table of other T values. T indexing is
%                 based on the upper index of T. For instance, T(i) would
%                 be equivalent to {super(i-1)}{T}{sub(i)}.

% Outputs:
%   DataStruct  - structure which represents the inverse homogeneous 
%                 transformation matrix.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last Edited 6/23/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Express a, alpha, d, and theta in terms of indexing in D
% a = D(index, 1);
% alpha = D(index, 2);
% d = D(index, 3);
% theta = D(index, 4);

% for inverse of homog matrix, [RT, -d; 0, 1] instead of [R, d; 0, 1];

DataStruct = HomogeneousTransform(index, D);

DataStruct = inv(DataStruct);

% Express inverse homogeneous transform matrix
% DataStruct = [cos(theta), sin(theta)*cos(alpha), sin(theta)*sin(alpha), -a; ...
%     -sin(theta), cos(theta)*cos(alpha), cos(theta)*sin(alpha), d*sin(alpha); ...
%     0, -sin(alpha), cos(alpha), -d*cos(alpha); ...
%     0, 0, 0, 1];

end