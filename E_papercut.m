% Cut graphing - Fold Type C
% Last edited 6/15/2021 by Lucien Peach

function [dataFoldE, m, lmax] = E_papercut(r, n, nl, ls, l1, l2, h0, h1, h2, alpha, beta)

% Counter used for data structure indexing
count = 1;

% Identify colors
orange = [1, 0.41, 0];
blue = [0, 0, 1];
black = [0, 0, 0];
red = [1, 0, 0];

% Begin by defining boundaries of sheet 
% -------------------------------------------------------------------

% This sheet boundary is irregular in that a small section on the right
% hand side is left blank, as it will be populated by the vertical fold
% region

% Overall boundary
boundarybottom = [0, 0; n*ls, 0];
boundaryleft = [0, 0; 0, h1 + 2*l2 + nl*2*l1 + h2 + l2];
boundarytop = [0, h1 + 2*l2 + nl*2*l1 + h2 + l2; n*ls, h1 + 2*l2 + nl*2*l1 + h2 + l2];
boundaryright = [n*ls, h1 + 2*l2 + nl*2*l1 + h2 + l2; n*ls, 0];
    
%[(n+1)*ls, h1 + 2*l2; (n+1)*ls, 0; 0, 0; 0, h1 + 2*l2 + nl*2*l1 + h2; ... 
%     (n+1)*ls, h1 + 2*l2 + nl*2*l1 + h2; (n+1)*ls, h1 + 2*l2 + nl*2*l1];

% Log data to structure and add to plot
dataFoldE(count).x = boundarybottom(:, 1);
dataFoldE(count).y = boundarybottom(:, 2);
dataFoldE(count).color = black;

figure()
plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
    dataFoldE(count).color)
hold on
set(gcf, 'color', 'w')

% Increase count
count = count + 1;

% Log data to structure and add to plot
dataFoldE(count).x = boundaryleft(:, 1);
dataFoldE(count).y = boundaryleft(:, 2);
dataFoldE(count).color = black;

plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
    dataFoldE(count).color)

% Increase count
count = count + 1;

% Log data to structure and add to plot
dataFoldE(count).x = boundarytop(:, 1);
dataFoldE(count).y = boundarytop(:, 2);
dataFoldE(count).color = black;

plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
    dataFoldE(count).color)

% Increase count
count = count + 1;

% Log data to structure and add to plot
dataFoldE(count).x = boundaryright(:, 1);
dataFoldE(count).y = boundaryright(:, 2);
dataFoldE(count).color = black;

plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
    dataFoldE(count).color)

% Increase count
count = count + 1;

% Log data to structure and add to plot
dataFoldE(count).x = boundaryleft(:, 1);
dataFoldE(count).y = boundaryleft(:, 2);
dataFoldE(count).color = blue;

plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
    dataFoldE(count).color)

% Increase count
count = count + 1;

% Log data to structure and add to plot
dataFoldE(count).x = boundaryright(:, 1);
dataFoldE(count).y = boundaryright(:, 2);
dataFoldE(count).color = blue;

plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
    dataFoldE(count).color)

% Increase count
% count = count + 1;

% Specify proximal line
% ------------------------------------------------------------------

% Store proximal line coordinates to array
% proximalx = [0; n*ls];
% proximaly = [h1; h1];
% dataFoldE(count).x = proximalx;
% dataFoldE(count).y = proximaly;
% dataFoldE(count).color = red;

% Plot proximal line
% plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
%     dataFoldE(count).color);

% The distal line will be overwritten by the top of the folded section, so
% we can ignore it for this purpose. 


% Specify boundaries of folded section
% ------------------------------------------------------------------

% Increase counter
count = count + 1;

% Store top horizontal fold lines to array
foldtopx = [0; n*ls];
foldtopy = [h1 + 2*l2 + nl*2*l1; h1 + 2*l2 + nl*2*l1];
dataFoldE(count).x = foldtopx;
dataFoldE(count).y = foldtopy;
dataFoldE(count).color = orange;

% Plot top boundary
plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
    dataFoldE(count).color)

% Increase counter
count = count + 1;

% Store bottom horizontal fold lines to array
foldbottomx = [0; n*ls];
foldbottomy = [h1 + 2*l2; h1 + 2*l2];
dataFoldE(count).x = foldbottomx;
dataFoldE(count).y = foldbottomy;
dataFoldE(count).color = orange;

% Plot bottom boundary of midsection
plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
    dataFoldE(count).color)

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
        orange_horizon(i-1, 2) = h1 + 2*l2 + (i/2)*2*l1;
        orange_horizon(i, 1) = n*ls;
        orange_horizon(i, 2) = h1 + 2*l2 + (i/2)*2*l1;
        
        % Store to array
        dataFoldE(count).x = orange_horizon(i-1:i, 1);
        dataFoldE(count).y = orange_horizon(i-1:i, 2);
        dataFoldE(count).color = orange;
        
        % Plot each line
        plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
            dataFoldE(count).color);
        
    end
end

% Specify horizontal line at l2:l2 boundary
% ------------------------------------------------------------------

% Increase counter
count = count + 1;

% Store x and y coordinates to array
lowerboundaryx = [0; n*ls];
lowerboundaryy = [h1 + l2; h1 + l2];
dataFoldE(count).x = lowerboundaryx;
dataFoldE(count).y = lowerboundaryy;
dataFoldE(count).color = blue;

% Begin plot. Plot line
plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
    dataFoldE(count).color);

% Specify vertical line pattern for non-folded region (blue lines)
% ------------------------------------------------------------------

% Bottom tube folds and graphing
bottomtube = zeros(2*(n-1), 2);

% Ignore side coordinates as these are graphed by the boundary section
for ii = 1:2:2*(n-1)
    
    % Increase count initially
    count = count + 1;
    
    % Indexing 
    index = ((ii-1)/2) + 1;
    
    % Populate array
    bottomtube(ii, 1) = index*ls;
    bottomtube(ii, 2) = 0;
    bottomtube(ii+1, 1) = index*ls;
    bottomtube(ii+1, 2) = h1 + 2*l2;
    
    % Log data to structure and add to plot. Plotting is sequential
    dataFoldE(count).x = bottomtube(ii:ii+1, 1);
    dataFoldE(count).y = bottomtube(ii:ii+1, 2);
    dataFoldE(count).color = blue;

    plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
        dataFoldE(count).color)
end

% Top tube folds and graphing
toptube = zeros(2*(n-1), 2);

% Ignore side folds as these are graphed by the boundary section
for jj = 1:2:2*(n-1)
    
    % Increase count initially
    count = count + 1;
    
    % Indexing 
    index = ((jj-1)/2) + 1;
    
    % Populate array
    toptube(jj, 1) = index*ls;
    toptube(jj, 2) = h1 + 2*l2 + nl*2*l1;
    toptube(jj+1, 1) = index*ls;
    toptube(jj+1, 2) = h1 + 2*l2 + nl*2*l1 + h2 + l2;
    
    % Log data to structure and add to plot. Plotting is sequential
    dataFoldE(count).x = toptube(jj:jj+1, 1);
    dataFoldE(count).y = toptube(jj:jj+1, 2);
    dataFoldE(count).color = blue;

    plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
        dataFoldE(count).color)
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
    midvert(kk-1, 2) = h1 + 2*l2;
    midvert(kk, 1) = index*ls;
    midvert(kk, 2) = h1 + 2*l2 + nl*2*l1;
    
    % Log data to structure and add to plot. Plotting is sequential
    dataFoldE(count).x = midvert(kk-1:kk, 1);
    dataFoldE(count).y = midvert(kk-1:kk, 2);
    dataFoldE(count).color = orange;

    plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
        dataFoldE(count).color)
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
            trianglefolds(aa+i, 2) = h1 + 2*l2 + (i-1)*l1;
 
        else % Odd
            
            % Determine x and y coordinates for point
            trianglefolds(aa+i, 1) = (index*ls);
            trianglefolds(aa+i, 2) = h1 + 2*l2 + (i-1)*l1;
         
        end
    end
    
    % Store to array
    dataFoldE(count).x = trianglefolds(aa+1:aa+(1 + 2*nl), 1);
    dataFoldE(count).y = trianglefolds(aa+1:aa+(1 + 2*nl), 2);
    dataFoldE(count).color = blue;
    
    % Plotting
    plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
        dataFoldE(count).color);
end

% Create horizontal segments (blue)
% ------------------------------------------------------------------

% Initialize no fold array
h_blue = zeros(2*(n)*nl, 2);

% Loop through to store [x,y] pairs for each line segment
for bb = 0:2*nl:2*nl*(n-1)
    
    % Increase count initially 
    count = count + 1;
    
    % Indexing 
    index = (bb/(2*nl)) + 1;
    
    % Add x and y coordinates for each segment
    for i = 1:2*nl
        
        % Coordinates will differ for "even" or "odd" points
        if rem(i, 2) == 0 % Even
            
            % Determine x and y coordinates for point
            h_blue(bb+i, 1) = (index*ls) - (l1 / tan(alpha));
            h_blue(bb+i, 2) = h1 + 2*l2 + (i-1)*l1;
 
        else % Odd
            
            % Determine x and y coordinates for point
            h_blue(bb+i, 1) = (index-1)*ls;
            h_blue(bb+i, 2) = h1 + 2*l2 + i*l1;
         
        end
    end
    
    % Storing to data array and plotting must take place in segments of two
    % points.
    for j = 0:2:(2*nl)-2
        dataFoldE(count).x = h_blue(bb+j+1:bb+j+2, 1);
        dataFoldE(count).y = h_blue(bb+j+1:bb+j+2, 2);
        dataFoldE(count).color = blue;
        
        plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
            dataFoldE(count).color)
    end 
end

% Create horizontal segments (orange)
% ------------------------------------------------------------------

% Initialize orange horizontal array
h_orange = zeros(2*(n)*nl, 2);

% Loop through to store [x,y] pairs for each line segment
for cc = 0:2*nl:2*nl*(n-1)
    
    % Increase count initially 
    count = count + 1;
    
    % Indexing 
    index = (cc/(2*nl)) + 1;
    
    % Add x and y coordinates for each segment
    for i = 1:2*nl
        
        % Coordinates will differ for "even" or "odd" points
        if rem(i, 2) == 0 % Even
            
            % Determine x and y coordinates for point
            h_orange(bb+i, 1) = index*ls;
            h_orange(bb+i, 2) = h1 + 2*l2 + (i-1)*l1;
 
        else % Odd
            
            % Determine x and y coordinates for point
            h_orange(bb+i, 1) = (index*ls) - (l1 / tan(alpha));
            h_orange(bb+i, 2) = h1 + 2*l2 + i*l1;
         
        end
    end
    
    % Storing to data array and plotting must take place in segments of two
    % points.
    for j = 0:2:(2*nl)-2
        dataFoldE(count).x = h_orange(bb+j+1:bb+j+2, 1);
        dataFoldE(count).y = h_orange(bb+j+1:bb+j+2, 2);
        dataFoldE(count).color = orange;
        
        plot(dataFoldE(count).x, dataFoldE(count).y, 'color', ...
            dataFoldE(count).color)
    end 
end

% Label the plot for clarity
title({
    ('Origami Schematic C for Provided Parameters:')
    ['[r = ' num2str(r) ', n = ' num2str(n) ', beta = ' num2str(beta) ...
    ', h0 = ' num2str(h0) ', nl = ' num2str(nl) ', dm = ' num2str(l2) ']']
    })

daspect([1 1 1])

m = 0;
lmax = h1 + 2*l2 + nl*2*l1 + h2 + l2;

close

end