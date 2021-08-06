% Finding the signed angle between two vector "a" and "b"  with respect to
% vector "a", around the normal vector "n". Outputs an angle "theta"
% Last edited 8/6/2021 by Wei-Hsi Chen

function [theta] = SignedAngle(a, b, n)

% Make sure all vectors are normalised
a = a/norm(a);
b = b/norm(b);
n = n/norm(n);

% find the signed angle
theta = atan2(dot(cross(a,b),n),dot(a,b));

