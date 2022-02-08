% Intersection Solver
% Last Edited 1/19/2021 by Lucien Peach

% Used to perform the mathematics necessary to find intersection point
% between a point and a line. Result is output as a position [x, y, z]
function [positionvector] = IntersectionSolver(plane, point, linevec)

% Begin by using the point and linevec information to create the values
% which will be used in parametric equations.
x0 = point(1);
y0 = point(2);
z0 = point(3);

a = linevec(1);
b = linevec(2);
c = linevec(3);

% Next, use plane to generate an equation for the plane. Recall that plane
% is a matrix.
normal = plane(:, 1);
center = plane(:, 2);

d = -1*(normal(1)*center(1) + normal(2)*center(2) + normal(3)*center(3));

% Thus, the equation of the plane is:
% normal(1)*x + normal(2)*y + normal(3)*z + d = 0;
% We can plug in our parametric equations for x, y, and z. This will allow
% us to determine the values of t, with which we can determine the
% intersection point.
const_sum = -(normal(1)*x0 + normal(2)*y0 + normal(3)*z0 + d);
t_sum = normal(1)*a + normal(2)*b + normal(3)*c;
t = const_sum / t_sum;

% Finally, determine intersection point
positionvector = zeros(3, 1);

positionvector(1) = x0 + a*t;
positionvector(2) = y0 + b*t;
positionvector(3) = z0 + c*t;

end
