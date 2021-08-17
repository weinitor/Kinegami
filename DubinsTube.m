% Dubins Tube Algorithm
% Last Updated 7/22/2021 by Lucien Peach
% Last Updated 7/27/2021 by Wei-Hsi Chen


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

tunit = t/norm(t);
infostruct(index).t = t;
infostruct(index).tunit = tunit;

% Determine values of wp and wd
wp = cross(ap, tunit);
wp = wp / norm(wp);

wd = cross(tunit, ad);
wd = wd / norm(wd);

% Define bm
[r_mat_p] = rotationalmatrix(wp, theta1);
bm = r_mat_p * bp;

% Determine values of phi1
% phi1 = atan2(norm(cross(bp, wp)), dot(bp, wp));
phi1 = SignedAngle(bp, wp, ap);

[r_mat_d] = rotationalmatrix(wd, theta2);

% Use results from rotational matrices, along with b_p and b_d, to
% determine the twist angle, alpha. 
bu = r_mat_d * bm;

% Use value of b_u to find alpha
% alpha = atan2(norm(cross(bu, bd)), dot(bu, bd));
alpha = SignedAngle(bu, bd, ad);

disp(alpha)

% Determine phi2
% phi2 = atan2(norm(cross(bm, wd)), dot(bm, wd)) - alpha;
phi2 = SignedAngle(bm, wd, tunit) - alpha;

% Elbow Fitting
[lengths, ls] = Origami_Elbow_creasedesign(r, n, phi1, theta1);

infostruct(index).ls = ls;
infostruct(index).r = r;

[dataFoldA, m, lmax] = Origami_Elbow_papercut(lengths, ls, n, infostruct(index).h1, ...
    infostruct(index).h2, r, phi1, theta1, mirror);

infostruct(index).m = m;
infostruct(index).lmax = lmax;
infostruct(index).n = n;
infostruct(index).type = dataFoldA;
infostruct(index).name = "Elbow";
infostruct(index).theta = theta1;
infostruct(index).dw = r*abs(tan(theta1 / 2));
infostruct(index).oc = Op(:, 4) + r*abs(tan(theta1 / 2))*ap;

% Twist Fitting
Tnorm = norm(t);

[x, l, ls] = Origami_Twist_creasedesign(r, n, 0.2*Tnorm, alpha);

infostruct(index+1).ls = ls;
infostruct(index+1).r = r;

[dataFoldB, m, lmax] = Origami_Twist_papercut(x, l, ls, n, infostruct(index+1).h1, ...
    infostruct(index+1).h2, r, 0.2*Tnorm, alpha);

infostruct(index+1).m = m;
infostruct(index+1).lmax = lmax;
infostruct(index+1).n = n;
infostruct(index+1).type = dataFoldB;
infostruct(index+1).name = "Twist";
infostruct(index+1).alpha = alpha;
infostruct(index+1).h = 0.2*Tnorm;

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
infostruct(index+2).h = 0.8*Tnorm;

% Elbow Fitting pt. 2
[lengths, ls] = Origami_Elbow_creasedesign(r, n, phi2, theta2);

infostruct(index+3).ls = ls;
infostruct(index+3).r = r;

[dataFoldD, m, lmax] = Origami_Elbow_papercut(lengths, ls, n, infostruct(index+3).h1, ...
    infostruct(index+3).h2, r, phi2, theta2, mirror);

infostruct(index+3).m = m;
infostruct(index+3).lmax = lmax;
infostruct(index+3).n = n;
infostruct(index+3).type = dataFoldD;
infostruct(index+3).name = "Elbow";
infostruct(index+3).theta = theta2;
infostruct(index+3).dw = r*abs(tan(theta2 / 2));
infostruct(index+3).oc = Op(:, 4) + r*abs(tan(theta1 / 2))*ap + t' + ...
    r*(abs(tan(theta1 / 2))+abs(tan(theta2 / 2)))*tunit';

% If the height of any segment is 0, edit so that lines are not printed
for i = index:index+3
    
    if infostruct(i).lmax == 0
        
        % Replaces x and y values with null, duplicate lines "erased"
        for j = 1:size(infostruct(i).type, 2)
            
            infostruct(i).type(j).x = [];
            infostruct(i).type(j).y = [];
            
        end
        
    end  
    
end

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
%         infostruct, i, msum, lmax_sum, 'triple');
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
% 
% close all


end
