% Dubins Path 3D Modeling
% Lasted edited 6/25/2021 by Dr. Cynthia Sung and Lucien Peach
% Lasted edited 7/27/2021 by Wei-Hsi Chen

% Outputs row vector and two radian values
% Make sure that Od and Op are normalized prior to entry
function [tmin, theta1min, theta2min] = solveDubins3d(r, Od, Op)

% input params
od = Od(:, 4).';
ad = Od(:, 1).';

op = Op(:, 4).';
ap = Op(:, 1).';

% normalize 'velocity vectors'
ad = ad / norm(ad);
ap = ap / norm(ap);

% solve for optimal t = [tx, ty, tz]
T0 = od-op+randn(size(od));

% Initialize vectors for error and solutions
Tsol = zeros(4, 3);
Terror = zeros(4, 3);

% provides individual solutions 
Tsol(1, :) = fsolve(@dubinspath1, T0); % + ... +
Tsol(2, :) = fsolve(@dubinspath2, T0); % - ... -
Tsol(3, :) = fsolve(@dubinspath3, T0); % + ... -
Tsol(4, :) = fsolve(@dubinspath4, T0); % - ... +

thetamatrix = zeros(4, 2);
ds = zeros(4, 1);
T_hat = zeros(4, 3);

for i = 1:4
    
    T_normalized = norm(Tsol(i, :));
    T_hat(i, :) = Tsol(i, :)/T_normalized;
    
% Find angle using arccos, the range of the output angle is (0, pi)   
%     thetamatrix(i, 1) = acos(dot(ap, T_hat(i, :)));
%     thetamatrix(i, 2) = acos(dot(ad, T_hat(i, :)));

% Find angle using atan2, the range of the output angle is (-pi, pi)   
    thetamatrix(i, 1) = atan2(norm(cross(ap, T_hat(i, :))),dot(ap, T_hat(i, :)));
    thetamatrix(i, 2) = atan2(norm(cross(T_hat(i, :), ad)),dot(T_hat(i, :), ad));
    
    ds(i) = r*thetamatrix(i, 1) + T_normalized + r*thetamatrix(i, 2);
    
end

[minds, indexds] = min(ds);

disp(minds);

% Output minimum Tunittor as well as minimum theta1 and theta2 values
tmin = Tsol(indexds, :);
theta1min = thetamatrix(indexds, 1);
theta2min = thetamatrix(indexds, 2);

% double check the error is low
Terror(1, :) = dubinspath1(Tsol(1, :));
Terror(2, :) = dubinspath2(Tsol(2, :));
Terror(3, :) = dubinspath3(Tsol(3, :));
Terror(4, :) = dubinspath4(Tsol(4, :));

% Display solutions and error
disp(Tsol)
disp(Terror)

% Plotting
% --------------------------------------------------------------------

% decompose path in relevant points for plotting
[Cd1, Pd1, Cp1, Pp1] = decomposePath1(Tsol(1,:));

% plotting
hold on
subplot(2, 2, 1)
plotCirclePath(op, od, Cp1, Pp1, Cd1, Pd1, ap, ad, 'b', 'r');

axis equal

% decompose path in relevant points for plotting
[Cd2, Pd2, Cp2, Pp2] = decomposePath2(Tsol(2,:));

% plotting
hold on
subplot(2, 2, 2)
plotCirclePath(op, od, Cp2, Pp2, Cd2, Pd2, ap, ad, 'b', 'r')

axis equal

% decompose path in relevant points for plotting
[Cd3, Pd3, Cp3, Pp3] = decomposePath3(Tsol(3,:));

% plotting
hold on
subplot(2, 2, 3)
plotCirclePath(op, od, Cp3, Pp3, Cd3, Pd3, ap, ad, 'b', 'r')

axis equal

% decompose path in relevant points for plotting
[Cd4, Pd4, Cp4, Pp4] = decomposePath4(Tsol(4,:));

% plotting
hold on
subplot(2, 2, 4)
plotCirclePath(op, od, Cp4, Pp4, Cd4, Pd4, ap, ad, 'b', 'r')

axis equal

% close all

% Functions
% ---------------------------------------------------------------------

    function [Cd, Pd, Cp, Pp] = decomposePath1(T) % + ... +
        % find relevant points based on T vector
        Tnorm = norm(T);
        Tunit = T/Tnorm;
        
        wp = cross(ap, cross(Tunit, ap));
        yp = cross(Tunit, cross(Tunit, ap));
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
        
        wd = -cross(ad, cross(Tunit, ad));
        yd = -cross(Tunit, cross(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

    function [Cd, Pd, Cp, Pp] = decomposePath2(T) % - ... -
        % find relevant points based on T vector
        Tnorm = norm(T);
        Tunit = T/Tnorm;
        
        wp = -cross(ap, cross(Tunit, ap));
        yp = -cross(Tunit, cross(Tunit, ap));
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
        
        wd = cross(ad, cross(Tunit, ad));
        yd = cross(Tunit, cross(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

    function [Cd, Pd, Cp, Pp] = decomposePath3(T) % + ... -
        % find relevant points based on T vector
        Tnorm = norm(T);
        Tunit = T/Tnorm;
        
        wp = cross(ap, cross(Tunit, ap));
        yp = cross(Tunit, cross(Tunit, ap));
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
        
        wd = cross(ad, cross(Tunit, ad));
        yd = cross(Tunit, cross(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

    function [Cd, Pd, Cp, Pp] = decomposePath4(T) % - ... +
        % find relevant points based on T vector
        Tnorm = norm(T);
        Tunit = T/Tnorm;
        
        wp = -cross(ap, cross(Tunit, ap));
        yp = -cross(Tunit, cross(Tunit, ap));
        Cp = op + r * wp / norm(wp); % center of starting circle
        Pp = Cp - r * yp / norm(yp); % point of leaving starting circle
        
        wd = -cross(ad, cross(Tunit, ad));
        yd = -cross(Tunit, cross(Tunit, ad));
        Cd = od + r * wd / norm(wd); % center of ending circle
        Pd = Cd - r * yd / norm(yd); % point of entering ending circle
    end

% 
% -------------------------------------------------------------------

%     function err = dubinspath(T)
%         % alternate formulation based just on vector math
%         [Cd, Pd, Cp, Pp] = decomposePath(T);
%         err = T - (Pp-Pd);
%     end

    function err = dubinspath1(T) % + ... +
        % equations to solve
        Tnorm = norm(T);
        Tunit = T/Tnorm;

        theta1 = atan2(norm(cross(ap,Tunit)),dot(ap,Tunit));
        theta2 = atan2(norm(cross(Tunit,ad)),dot(Tunit,ad));

        err = T + ...
            r*(tan(theta1/2) + tan(theta2/2))*Tunit + ...
            r*(tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
            - (od-op);
    end

    function err = dubinspath2(T) % - ... -
        % equations to solve
        Tnorm = norm(T);
        Tunit = T/Tnorm;

        theta1 = atan2(norm(cross(ap,Tunit)),dot(ap,Tunit));
        theta2 = atan2(norm(cross(Tunit,ad)),dot(Tunit,ad));

        err = T - ...
            r*(tan(theta1/2) + tan(theta2/2))*Tunit - ...
            r*(tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
            - (od-op);
    end

    function err = dubinspath3(T) % + ... -
        % equations to solve
        Tnorm = norm(T);
        Tunit = T/Tnorm;

        theta1 = atan2(norm(cross(ap,Tunit)),dot(ap,Tunit));
        theta2 = atan2(norm(cross(Tunit,ad)),dot(Tunit,ad));

        err = T + ...
            r*(tan(theta1/2) - tan(theta2/2))*Tunit + ...
            r*(tan(theta1/2)*ap - tan(theta2/2)*ad) + ...
            - (od-op);
    end

    function err = dubinspath4(T) % - ... +
        % equations to solve
        Tnorm = norm(T);
        Tunit = T/Tnorm;

        theta1 = atan2(norm(cross(ap,Tunit)),dot(ap,Tunit));
        theta2 = atan2(norm(cross(Tunit,ad)),dot(Tunit,ad));

        err = T + ...
            r*(-tan(theta1/2) + tan(theta2/2))*Tunit + ...
            r*(-tan(theta1/2)*ap + tan(theta2/2)*ad) + ...
            - (od-op);
    end

    % Plotting without overriding
    function plotCirclePath(Op, Od, Cp, Pp, Cd, Pd, ap, ad, col1, col2)
        
        hold on
        plot3(Op(1),Op(2),Op(3), ['*' col2])
        plot3(Op(1)+[0,r*ap(1)],Op(2)+[0,r*ap(2)],Op(3)+[0,r*ap(3)],['-' col2])
        plot3(Cp(1),Cp(2),Cp(3), ['*' col2])
        plot3(Pp(1),Pp(2),Pp(3),'*g')
        
        n = cross(Op-Cp, Pp-Cp);
        v = null(n);
        theta = linspace(0,2*pi,50);
        circle_pts = bsxfun(@plus, Cp(:), r * (v(:,1)*cos(theta) + v(:,2)*sin(theta)));
        plot3(circle_pts(1,:),circle_pts(2,:),circle_pts(3,:), ['--' col2]);
        
        plot3(Od(1),Od(2),Od(3), ['*' col1])
        plot3(Od(1)+[0,r*ad(1)],Od(2)+[0,r*ad(2)],Od(3)+[0,r*ad(3)],['-' col1])
        plot3(Cd(1),Cd(2),Cd(3), ['*' col1])
        plot3(Pd(1),Pd(2),Pd(3),'*g')
        
        n = cross(Od-Cd, Pd-Cd);
        v = null(n);
        theta = linspace(0,2*pi,50);
        circle_pts = bsxfun(@plus, Cd(:), r * (v(:,1)*cos(theta) + v(:,2)*sin(theta)));
        plot3(circle_pts(1,:),circle_pts(2,:),circle_pts(3,:), ['--' col1]);
        
        
        plot3([Pd(1) Pp(1)], [Pd(2) Pp(2)], [Pd(3) Pp(3)],'-g')
    end

end