% No-Intersect Algorithm
% Last Edited 8/24/2021 by Lucien Peach

function [TransformStruct] = Intersection(TransformStruct, index, r)

% ----------- Old Code - Keep in Case of Intersection Checks ------------

% % Check if straight line path intersects previous spheres:
% 
% % Define variables that will be used in calculations
% syms t
%
% % Determine directionality vectors
% vecx = CoordVec(1); % (x2 - x1)
% vecy = CoordVec(2); % (y2 - y1)
% vecz = CoordVec(3); % (z2 - z1)
%
% 
% 
% % Define x1, x2; y1, y2; z1, z2
% x1 = TransformStruct(index).Oc(1, 4);
% x2 = TransformStruct(index+1).Oc(1, 4);
% y1 = TransformStruct(index).Oc(2, 4);
% y2 = TransformStruct(index+1).Oc(2, 4);
% z1 = TransformStruct(index).Oc(3, 4);
% z2 = TransformStruct(index).Oc(3, 4);
% 
% % Determine a, b, and c values from quadratic equation
% a = (vecx)^2 + (vecy)^2 + (vecz)^2;
% 
% % Use initial point and directionality vectors to determine x, y, and z
% % values, which will be plugged into the sphere equation
% % x = TransformStruct(index).Oc(1, 4) + vecx*t;
% % y = TransformStruct(index).Oc(2, 4) + vecy*t;
% % z = TransformStruct(index).Oc(3, 4) + vecz*t;
% 
% % Determine a, b, and c values from quadratic equation
% sphere_eqn = ((x - TransformStruct(index).oinew(1))^2) + ((y - ...
%     TransformStruct(index).oinew(2))^2) + ((z - ...
%     TransformStruct(index).oinew(3))^2) - 16 == 0;
% 
% % Determine coordinates 

% --------------------------- Optimization ------------------------------

% Find minimal distance between line and boundary sphere center using
% parametrized line equations

% Define z-axis, along which point can move
zzz = TransformStruct(index).Oc(:, 3);

% Define total goal dist
dist = 4*r + TransformStruct(index+2).rnew;

% To find parametric equation of line, take center points from consecutive
% spheres
A = TransformStruct(index+2).oilarge;
B = TransformStruct(index+1).Oc(:, 4);
C = TransformStruct(index).Oc(:, 4);

CoordVec = C - B;

% Determine direction vector
d = CoordVec / norm(CoordVec);

% Determine vector from previous bounding sphere center and any point on
% the line (we can use the point of sphere (index+1))
v = A - B;

% Dot Product
t = v*d;

% Projection of point A onto line and final distance calculation
P = B + (t*d);
boundingdist = norm(P - A);

% See if minimum distance criteria is met or not
if boundingdist < 4*r
    
    % Run optimization
    % Statement for minimization
    obj = @(delta) norm(C + zzz*delta - B);
    
    % Initial Guess
    delta0 = 1;
    
    % Inequality requirements:
    % norm(B
    
    % Nonlinearity constraint
    function [c, ceq] = nlcon(delta)
        c = dist - norm()
        ceq = [];
    
    end
    
end

    % Inequality requirements:
    % norm(DataStruct(index-1).oi + zzz*delta - DataStruct(index).oi) >= dist
    % Nonlinear

        % Nonlinearity constraint
        function [c, ceq] = nlcon(delta)
            c = dist - norm(DataStruct(index-1).oi + zzz*delta - DataStruct(index).oilarge);
            ceq = [];
        end

    % Express other parameters as empty cells (not used)
    A = [];
    B = [];
    Aeq = [];
    Beq = [];
    LB = [];
    UB = [];
    nonlincon = @nlcon;

    % Output optimal value for delta
    delta_opt = fmincon(obj, delta0, A, B, Aeq, Beq, LB, UB, nonlincon);

end


% 2) Determine equation of sphere
% 3) Determine equation of line between sphere centers
% 4) Test intersection of sphere equation and line equation
% 5) If intersection, proceed. If no intersection, exit function


% Vary delta such that:
% 1) Shortest distance between line is 4r
% 2) distance between two spheres is minimal 


% If so, move the next sphere along its z axis until it reaches a spot that
% is at least 4r away from the current sphere boundaries but also creates a
% straight path to the previous sphere that does not intersect 

% We want the shortest possible distance displaced

