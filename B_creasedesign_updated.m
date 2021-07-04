% Basic Crease Design - Fold Type B
% Last edited 6/4/2021 by Lucien Peach
% Modified by Wei-Hsi 06/22/2021

% Function declaration
function [x, l, ls] = B_creasedesign_updated(r, n, h, alpha)

% Specify delta value for midsection height calculation
delta = pi*((n-2)/(2*n));

% Side length
ls = 2*r*cos(delta);

% Specify parameter for 1/2 base of midsection triangles
% x = r * mod(alpha, (2*pi)/n);
a = mod(alpha, (2*pi)/n);
x = ls/2 * (1-cos(a)+cot(pi/n)*sin(a));

% Length of midsection
l = sqrt((h^2) + (ls*csc(pi/n)*sin(pi/n-a/2)*sin(a/2))^2);

% this is how to calculate the marker
% m = floor(alpha/((2*pi)/n))*ls + x;

end
