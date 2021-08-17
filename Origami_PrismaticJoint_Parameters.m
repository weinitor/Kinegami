% Crease pattern parameters - Origami Prismatic joint
% Last edited 6/15/2021 by Lucien Peach

% Function declaration
function [ls, l1, h0, dm, alpha] = Origami_PrismaticJoint_Parameters(r, n, beta, d0, nl)

% Specify delta value for side length calculation
delta = pi*((n-2)/(2*n));

% Side length
ls = 2*r*cos(delta);

% Determine value of alpha for use in angular measurements
alpha = (pi - (2*pi*cos(beta)/n))/2;

% Determine h0
h0 = d0/nl;

% Determine value of l1, the height of the midsection 
l1 = h0*csc(beta)/2;

% Determine max length quantity
dm = nl * 2* l1;

% Display an error prompt if l2 is too small. Otherwise, return a value
% that can be checked to allow the code to proceed 
% if l2 > req
%     alert = 0;
% else
%     alert = msgbox('Choose a Larger Value for l2', 'Error', 'error');
% end

end
