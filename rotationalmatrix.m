% Rotation Matrix
% Last edited 6/16/2021 by Lucien Peach

function [r_mat] = rotationalmatrix(w, theta)

wx = w(1);
wy = w(2);
wz = w(3);

c = cos(theta);
s = sin(theta);

% Initialize r_mat
r_mat = zeros(3, 3);

% Populate r_mat, the 3x3 rotational matrix

% Column 1
r_mat(1, 1) = (wx^2)*(1 - c) + c;
r_mat(2, 1) = wx*wy*(1 - c) + (wz*s);
r_mat(3, 1) = wx*wz*(1 - c) - (wy*s);

% Column 2
r_mat(1, 2) = wx*wy*(1 - c) - (wz*s);
r_mat(2, 2) = (wy^2)*(1 - c) + c;
r_mat(3, 2) = wy*wz*(1 - c) + (wx*s);

% Column 3
r_mat(1, 3) = wx*wz*(1 - c) + (wy*s);
r_mat(2, 3) = wy*wz*(1 - c) - (wx*s);
r_mat(3, 3) = (wz^2)*(1 - c) + c;

end