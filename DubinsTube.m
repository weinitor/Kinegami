% Dubins Tube Algorithm
% Last Updated 6/27/2021 by Lucien Peach

function [infostruct] = DubinsTube(r, n, Op, Od, infostruct, index, mirror)

% Op and Od are both 3x4 matrices

% We know that Ou = {au, bu, cu, ou} - 3x4 matrix of vec/norm(vec)

% Thus,
ap = Op(:, 1);
ap = ap / norm(ap);

ad = Od(:, 1);
ad = ad / norm(ad);

bp = Op(:, 2);
bp = bp / norm(bp);

bd = Od(:, 2);
bd = bd / norm(bd);

% Run solution solver to compare four sets of solutions and find t with the
% shortest path. output this t and the corresponding theta1 and theta2
[t, theta1, theta2] = solveDubins3d(r, Od, Op);

tvec = t/norm(t);

% Determine values of wp and wd
wp = cross(ap, tvec);
% wp = wp / norm(wp);
% wd = cross(ad, tvec);
wd = cross(tvec, ad);
% wd = wd / norm(wd);

% Define bm
[r_mat_p] = rotationalmatrix(wp, theta1);
bm = r_mat_p * bp;

% Determine values of phi1
phi1 = atan2(norm(cross(bp, wp)), dot(bp, wp));

[r_mat_d] = rotationalmatrix(wd, theta2);

% Use results from rotational matrices, along with b_p and b_d, to
% determine the twist angle, alpha. 
bu = r_mat_d * bm;

% Use value of b_u to find alpha
alpha = atan2(norm(cross(bu, bd)), dot(bu, bd));

disp(alpha)

% Determine phi2
phi2 = atan2(norm(cross(bm, wd)), dot(bm, wd)) - alpha;

% Joint Definitions and Fittings
h1 = 0;
h2 = 0;

% Elbow Fitting
[lengths, ls] = A_creasedesign(r, n, phi1, theta1);

infostruct(index).ls = ls;
infostruct(index).r = r;

[dataFoldA, m, lmax] = A_papercut(lengths, ls, n, h1, h2, r, phi1, theta1, mirror);

infostruct(index).m = m;
infostruct(index).lmax = lmax;
infostruct(index).n = n;
infostruct(index).type = dataFoldA;
infostruct(index).name = "Elbow";

% Twist Fitting
Tnorm = norm(t);

[x, l, ls] = B_creasedesign_updated(r, n, 0.2*Tnorm, alpha);

infostruct(index+1).ls = ls;
infostruct(index+1).r = r;

[dataFoldB, m, lmax] = B_papercut(x, l, ls, n, h1, h2, r, 0.2*Tnorm, alpha);

infostruct(index+1).m = m;
infostruct(index+1).lmax = lmax;
infostruct(index+1).n = n;
infostruct(index+1).type = dataFoldB;
infostruct(index+1).name = "Twist";

% Tube
[ls] = Default_creasedesign(r, n);

infostruct(index+2).ls = ls;
infostruct(index+2).r = r;

% Outputs graphing for default tube
[dataFoldTube, m, lmax] = Default_papercut(n, ls, 0.8*Tnorm, r);

infostruct(index+2).m = m;
infostruct(index+2).lmax = lmax;
infostruct(index+2).n = n;
infostruct(index+2).type = dataFoldTube;
infostruct(index+2).name = "Tube";

% Elbow Fitting pt. 2
[lengths, ls] = A_creasedesign(r, n, phi2, theta2);

infostruct(index+3).ls = ls;
infostruct(index+3).r = r;

[dataFoldD, m, lmax] = A_papercut(lengths, ls, n, h1, h2, r, phi2, theta2, mirror);

infostruct(index+3).m = m;
infostruct(index+3).lmax = lmax;
infostruct(index+3).n = n;
infostruct(index+3).type = dataFoldD;
infostruct(index+3).name = "Elbow";

% % Add field for tracking lmaxnet
% infostruct(1).lmaxnet = infostruct(1).lmax;
% 
% msum = 0;
% lmax_sum = 0;
% 
% % figure()
% set(gcf, 'color', 'w')
% hold on
% 
% % Loop through indices to plot
% for i = index:index+3
%     
%     [msum, lmax_sum] = DataFoldAppend(infostruct(i).type, ...
%     infostruct, i, msum, lmax_sum);
%     
% end
% 
% msum = 0;
% lmax_sum = 0;
% 
% % Loop through indices to plot duplication
% for i = index:index+3
%     
%     [msum, lmax_sum] = DataFoldDuplicate(infostruct(i).type, ...
%         infostruct, i, msum, lmax_sum);
%     
% end
% 
% for i = 1:size(infostruct, 2)
% 
%     if i == 1
%         lmaxtotal = infostruct(i).lmax;
%     else
%         lmaxtotal = infostruct(i).lmax + lmaxtotal;
%     end
% 
% end
% 
% % Figure adjustments
% daspect([1, 1, 1])
% axis off

% close all


end