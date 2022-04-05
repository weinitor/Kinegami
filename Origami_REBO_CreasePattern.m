function [dataFoldREBO, m, lmax] = Origami_REBO_CreasePattern(n, nl, ls, l1, h1, h2, alpha)
% ORIGAMI_REBO_CREASEPATTERN - Generates a crease pattern for the origami
% REBO joint.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last edited 4/3/2022
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.txt for detail.


% Counter used for data structure indexing
count = 1;

% Identify colors
orange = [1, 0.41, 0];
blue = [0, 0, 1];
black = [0, 0, 0];

% Begin by defining boundaries of sheet 
% -------------------------------------------------------------------

% We will be accommodating for a fold-over panel as well. 

% Overall boundary
boundarybottom = [0, 0; n*ls, 0];
boundaryleft = [0, 0; 0, h1 + nl*2*l1 + h2];
boundarytop = [0, h1 + nl*2*l1 + h2; n*ls, h1 + nl*2*l1 + h2];
boundaryright = [n*ls, h1 + nl*2*l1 + h2; n*ls, 0];

% Log data to structure and add to plot
dataFoldREBO(count).x = boundarybottom(:, 1);
dataFoldREBO(count).y = boundarybottom(:, 2);
dataFoldREBO(count).color = black;

% Increase count
count = count + 1;

% Log data to structure and add to plot
dataFoldREBO(count).x = boundaryleft(:, 1);
dataFoldREBO(count).y = boundaryleft(:, 2);
dataFoldREBO(count).color = black;

% Increase count
count = count + 1;

% Log data to structure and add to plot
dataFoldREBO(count).x = boundarytop(:, 1);
dataFoldREBO(count).y = boundarytop(:, 2);
dataFoldREBO(count).color = black;

% Increase count
count = count + 1;

% Log data to structure and add to plot
dataFoldREBO(count).x = boundaryright(:, 1);
dataFoldREBO(count).y = boundaryright(:, 2);
dataFoldREBO(count).color = black;

% Specify intermediary horizontal lines for nl > 1
% ------------------------------------------------------------------
if nl > 1
    
    % Initialization of horizontal lines vector
    orange_horizon = zeros(2*(nl-1), 2);
    
    for i = 2:2:(nl*2)-2
        
        % Increase count initially
        count = count + 1;
        
        % Populate x and y values for each horizontal line
        orange_horizon(i-1, 1) = 0;
        orange_horizon(i-1, 2) = h1 + (i/2)*2*l1;
        orange_horizon(i, 1) = n*ls;
        orange_horizon(i, 2) = h1 + (i/2)*2*l1;
        
        % Store to array
        dataFoldREBO(count).x = orange_horizon(i-1:i, 1);
        dataFoldREBO(count).y = orange_horizon(i-1:i, 2);
        dataFoldREBO(count).color = orange;
        
    end
end

% Specify vertical line pattern (orange)
% ------------------------------------------------------------------

% Bottom tube folds and graphing
midvert = zeros(2*(n), 2);

% Ignore side folds as these are graphed by the boundary section
for kk = 2:2:2*(n)
    
    % Increase count initially
    count = count + 1;
    
    % Indexing 
    index = ((kk-2)/2) + 1;
    
    % Populate array
    midvert(kk-1, 1) = index*ls;
    midvert(kk-1, 2) = h1;
    midvert(kk, 1) = index*ls;
    midvert(kk, 2) = h1 + nl*2*l1;
    
    % Log data to structure and add to plot. Plotting is sequential
    dataFoldREBO(count).x = midvert(kk-1:kk, 1);
    dataFoldREBO(count).y = midvert(kk-1:kk, 2);
    dataFoldREBO(count).color = orange;

end

% Create triangular (blue)
% ------------------------------------------------------------------

% Initialize array for triangular folds
ridge = (1 + 2*nl)*(n);
trianglefolds = zeros(ridge, 2);

% Specify x and y coordinates for each column of triangles
for aa = 0:(1 + 2*nl):ridge - (1 + 2*nl)
    
    % Increase count initially
    count = count + 1;
    
    % Indexing
    index = 1 + (aa/(1 + 2*nl)); 
    
    % Create a variable that will reset each loop. These will be used to
    % determine where in the column each "even" or "odd" point resides.
    % tower_even = 0;
    % tower_odd = 0;
    
    % Add x and y coordinates for column of triangles
    for i = 1:(1 + 2*nl)
        
        % Coordinates will differ for "even" or "odd" points
        if rem(i, 2) == 0 % Even
            
            % Determine x and y coordinates for point
            trianglefolds(aa+i, 1) = (index*ls) - (l1 / tan(alpha));
            trianglefolds(aa+i, 2) = h1 + (i-1)*l1;
 
        else % Odd
            
            % Determine x and y coordinates for point
            trianglefolds(aa+i, 1) = (index*ls);
            trianglefolds(aa+i, 2) = h1 + (i-1)*l1;
         
        end
    end
    
    % Store to array
    dataFoldREBO(count).x = trianglefolds(aa+1:aa+(1 + 2*nl), 1);
    dataFoldREBO(count).y = trianglefolds(aa+1:aa+(1 + 2*nl), 2);
    dataFoldREBO(count).color = blue;
    
end

% Create horizontal segments (blue)
% ------------------------------------------------------------------

% Initialize no fold array
h_blue = zeros(2*(n)*nl, 2);

% Loop through to store [x,y] pairs for each line segment
for bb = 0:2*nl:2*nl*(n-1)
    
    % Indexing 
    index = (bb/(2*nl)) + 1;
    
    % Add x and y coordinates for each segment
    for i = 1:2*nl
        
        % Coordinates will differ for "even" or "odd" points
        if rem(i, 2) == 0 % Even
            
            % Determine x and y coordinates for point
            h_blue(bb+i, 1) = (index*ls) - (l1 / tan(alpha));
            h_blue(bb+i, 2) = h1 + (i-1)*l1;
 
        else % Odd
            
            % Determine x and y coordinates for point
            h_blue(bb+i, 1) = (index-1)*ls;
            h_blue(bb+i, 2) = h1 + i*l1;
         
        end
    end
    
    % Storing to data array and plotting must take place in segments of two
    % points.
    for j = 0:2:(2*nl)-2
        
        % Increase count initially 
        count = count + 1;
        
        dataFoldREBO(count).x = h_blue(bb+j+1:bb+j+2, 1);
        dataFoldREBO(count).y = h_blue(bb+j+1:bb+j+2, 2);
        dataFoldREBO(count).color = blue;
        
    end 
end

% Create horizontal segments (orange)
% ------------------------------------------------------------------

% Initialize orange horizontal array
h_orange = zeros(2*(n)*nl, 2);

% Loop through to store [x,y] pairs for each line segment
for cc = 0:2*nl:2*nl*(n-1)
    
    % Indexing 
    index = (cc/(2*nl)) + 1;
    
    % Add x and y coordinates for each segment
    for i = 1:2*nl
        
        % Coordinates will differ for "even" or "odd" points
        if rem(i, 2) == 0 % Even
            
            % Determine x and y coordinates for point
            h_orange(cc+i, 1) = index*ls;
            h_orange(cc+i, 2) = h1 + (i-1)*l1;
 
        else % Odd
            
            % Determine x and y coordinates for point
            h_orange(cc+i, 1) = (index*ls) - (l1 / tan(alpha));
            h_orange(cc+i, 2) = h1 + i*l1;
         
        end
    end
    
    % Storing to data array and plotting must take place in segments of two
    % points.
    for j = 0:2:(2*nl)-2
        
        % Increase count initially 
        count = count + 1;
        
        dataFoldREBO(count).x = h_orange(cc+j+1:cc+j+2, 1);
        dataFoldREBO(count).y = h_orange(cc+j+1:cc+j+2, 2);
        dataFoldREBO(count).color = orange;

    end 
end

m = 0;
lmax = h1 + nl*2*l1 + h2;

end