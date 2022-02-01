% Dubins Path Plotting
% Last edited 7/22/2021 by Lucien Peach
% Last edidt 7/27/2021 by Wei-Hsi Chen

function [TransformStruct] = DubinsPlot(TransformStruct, infostruct, index, i)

% TransformStruct provides information about frames
% infostruct provides information about fold parameters
% i used for TransformStruct indexing
% index used for infostruct indexing

% Using same index scheme as for DubinsTube: 
% index = elbow1
% index+1 = twist
% index+2 = tube
% index+3 = elbow2

tunit = infostruct(index).tunit;
t = infostruct(index).t;

% Origin point (Od)
TransformStruct(i).path(:, 1) = TransformStruct(i).Od(:, 4);

% Constructing the path via the parameters of the DubinsTube
% Elbow1 Parameters and Path Plotting
% ------------------------------------------------------------------------

% Pull ap from initial frame
ap = TransformStruct(i).Od(:, 1);
ap = ap/norm(ap);

% Determine ad based on theta - rotation matrix depends on axis of rotation
% First, we must find w (w must be a normal vector)
wp = cross(ap, tunit');
wp = wp/norm(wp);

% Use ap, RotationalMatrix, and w to determine ad for first elbow
ad = RotationalMatrix(wp, infostruct(index).theta)*ap;
ad = ad/norm(ad);

% Find the vector that points to the instance center of rotation for the
% fisrt elbow
o2c1 = -cross(ap, wp);
o2c1 = o2c1/norm(o2c1);
% center1 = TransformStruct(i).path(:, 1) + infostruct(index).r*o2c1;

% Factor in h1 and h2 offsets to determine positions for plotting
TransformStruct(i).path(:, 2) = (infostruct(index).h1 + infostruct(index).dw)*ap ...
    + TransformStruct(i).path(:, 1);
TransformStruct(i).path(:, 3) = (infostruct(index).h2 + infostruct(index).dw)*ad ...
    + TransformStruct(i).path(:, 2);


% Twist Parameters and Path Plotting
% ------------------------------------------------------------------------

% Frame direction remains consistent - ad from elbow1
% For purposes of initial iteration, twist functions similarly to tube but
% with an adjusted height calculation
TransformStruct(i).path(:, 4) = (infostruct(index+1).h1 + infostruct(index+1).h2 + ...
    infostruct(index+1).h)*ad + TransformStruct(i).path(:, 3);


% Tube Parameters and Path Plotting
% ------------------------------------------------------------------------

% Frame direction remains consistent - ad from elbow1
TransformStruct(i).path(:, 5) = (infostruct(index+2).h1 + infostruct(index+2).h2 + ...
    infostruct(index+2).h)*ad + TransformStruct(i).path(:, 4);


% Elbow2 Parameters and Path Plotting
% ------------------------------------------------------------------------

% Pull ap from elbow1 ad
ap2 = ad;
a_next = TransformStruct(i+1).Op(:, 1);
a_next = a_next/norm(a_next);

% Determine new ad based on theta, (w must be a normal vector)
om2od = TransformStruct(i+1).Op(:, 4) - TransformStruct(i).path(:, 5);
om2od = om2od/norm(om2od);
% wp2 = cross(ap2, a_next);
wp2 = cross(ap2, om2od);
wp2 = wp2/norm(wp2);
ad2 = RotationalMatrix(wp2, infostruct(index+3).theta) * ap2;
ad2 = ad2/norm(ad2);

% Find the vector that points to the instance center of rotation for the
% second elbow
o2c2 = -cross(ap2,wp2);
o2c2 = o2c2/norm(o2c2);
center2 = TransformStruct(i).path(:, 5) + infostruct(index).r*o2c1;

% Factor in h1 and h2 offsets
TransformStruct(i).path(:, 6) = (infostruct(index+3).h1 + infostruct(index+3).dw)*ap2 ...
    + TransformStruct(i).path(:, 5);
TransformStruct(i).path(:, 7) = (infostruct(index+3).h2 + infostruct(index+3).dw)*ad2 ...
    + TransformStruct(i).path(:, 6);

% Wei: What is this section exactly???
% Constructing the path via the parameters of the DubinsPath
% Newpath
% ------------------------------------------------------------------------
TransformStruct(i).newpath(:, 1) = TransformStruct(i).Od(:, 4);

% TransformStruct(i).newpath(:, 2) = (infostruct(index).h1 + infostruct(index).dw)*ap ...
%     + TransformStruct(i).newpath(:, 1);
TransformStruct(i).newpath(:, 2) = infostruct(index).oc;

% TransformStruct(i).newpath(:, 3) = TransformStruct(i).newpath(:, 2) + ...
%     t.'+ (infostruct(index+3).dw + infostruct(index).dw)*tunit.';
TransformStruct(i).newpath(:, 3) = infostruct(index+3).oc;

% TransformStruct(i).newpath(:, 4) = TransformStruct(i+1).Op(:, 4);
TransformStruct(i).newpath(:, 4) = TransformStruct(i).newpath(:, 3) + ...
    TransformStruct(i+1).Op(:, 1)*infostruct(index+3).dw;

% -------------------------
% Plot the two circular arcs
dAng1 = infostruct(index).theta/50;
circle_pts1 = [];
for j = 1:50
    circle_pts1 = [circle_pts1, TransformStruct(i).path(:, 1) + ...
        infostruct(index).r * (eye(3) - RotationalMatrix(wp, dAng1*j)) * o2c1];
end

if infostruct(index+3).theta < 0
    dAng2 = (2*pi + infostruct(index+3).theta)/50;
else
    dAng2 = infostruct(index+3).theta/50;
end
circle_pts2 = [];
for j = 0:50
    circle_pts2 = [circle_pts2, TransformStruct(i).path(:, 5) + ...
        infostruct(index+3).r * (eye(3) - RotationalMatrix(wp2, dAng2*j)) * o2c2];
end

CSCDubinspath = [circle_pts1, circle_pts2];

% Should now create a "path" represented by segments 
% handle = plot3(TransformStruct(i).path(1, 3:5), TransformStruct(i).path(2, 3:5), ...
%     TransformStruct(i).path(3, 3:5), 'LineWidth', 4, 'Color', 'k');
% handle.Annotation.LegendInformation.IconDisplayStyle = 'off';

% handle = plot3(circle_pts1(1,:),circle_pts1(2,:),circle_pts1(3,:),...
%     TransformStruct(i).path(1, 3:5), TransformStruct(i).path(2, 3:5), ...
%     TransformStruct(i).path(3, 3:5), ...
%     circle_pts2(1,:),circle_pts2(2,:),circle_pts2(3,:),...
%     'LineWidth', 4, 'Color', 'k');
% handle.Annotation.LegendInformation.IconDisplayStyle = 'off';

% handle2 = plot3(TransformStruct(i).newpath(1, 1:4), TransformStruct(i).newpath(2, 1:4), ...
%     TransformStruct(i).newpath(3, 1:4), 'LineWidth', 3, 'Color', 'k');
% handle2.Annotation.LegendInformation.IconDisplayStyle = 'off';

% handle2 = plot3(circle_pts1(1,:),circle_pts1(2,:),circle_pts1(3,:),...
%     circle_pts2(1,:),circle_pts2(2,:),circle_pts2(3,:),'LineWidth', 4, 'Color', 'k');

handle = plot3(CSCDubinspath(1,:), CSCDubinspath(2,:), CSCDubinspath(3,:),...
    'LineWidth', 4, 'Color', 'k');
handle.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
