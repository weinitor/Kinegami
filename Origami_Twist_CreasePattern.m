% Graph for crease pattern - Origami twist fitting
% Last edited 6/7/2021 by Lucien Peach

function [dataFoldB, m, lmax] = Origami_Twist_CreasePattern(x, l, ls, n, h1, h2, r, h, alpha)

% Counter used for data structure indexing
count = 1;

% Identify colors
orange = [1, 0.41, 0];
blue = [0, 0, 1];
black = [0, 0, 0];

% Overall boundary
boundarybottom = [0, 0; n*ls, 0];
boundaryleft = [0, 0; 0, h1 + l + h2];
boundarytop = [0, h1 + l + h2; n*ls, h1 + l + h2];
boundaryright = [n*ls, h1 + l + h2; n*ls, 0];

% Log data to structure and add to plot
dataFoldB(count).x = boundarybottom(:, 1);
dataFoldB(count).y = boundarybottom(:, 2);
dataFoldB(count).color = black;

% Increase counter
count = count + 1;

% Log data to structure and add to plot
dataFoldB(count).x = boundaryleft(:, 1);
dataFoldB(count).y = boundaryleft(:, 2);
dataFoldB(count).color = black;

% Increase counter
count = count + 1;

% Log data to structure and add to plot
dataFoldB(count).x = boundarytop(:, 1);
dataFoldB(count).y = boundarytop(:, 2);
dataFoldB(count).color = black;

% Increase counter
count = count + 1;

% Log data to structure and add to plot
dataFoldB(count).x = boundaryright(:, 1);
dataFoldB(count).y = boundaryright(:, 2);
dataFoldB(count).color = black;

% Log data to structure and add to plot
% dataFoldB(count).x = boundaryleft(:, 1);
% dataFoldB(count).y = boundaryleft(:, 2);
% dataFoldB(count).color = blue;
% 
% plot(dataFoldB(count).x, dataFoldB(count).y, 'color', ...
%     dataFoldB(count).color)
% 
% % Increase counter
% count = count + 1;
% 
% % Log data to structure and add to plot
% dataFoldB(count).x = boundaryright(:, 1);
% dataFoldB(count).y = boundaryright(:, 2);
% dataFoldB(count).color = blue;
% 
% plot(dataFoldB(count).x, dataFoldB(count).y, 'color', ...
%     dataFoldB(count).color)
% 
% count = count + 1;

% plot(dataPaper{1}(1).x, dataPaper{1}(1).y)

% Special case exceptions
if x > 10^-4 && x < ls - 10^-4

    % Increase counter
    count = count + 1;

    % Begin by specifying top and bottom parameters (orange lines)
    topsegmentx = [0; n*ls];
    topsegmenty = [h1+l; h1+l];
    dataFoldB(count).x = topsegmentx;
    dataFoldB(count).y = topsegmenty;
    dataFoldB(count).color = orange;

    % Increase counter
    count = count+1;

    bottomsegmentx = [0; n*ls];
    bottomsegmenty = [h1; h1];
    dataFoldB(count).x = bottomsegmentx;
    dataFoldB(count).y = bottomsegmenty;
    dataFoldB(count).color = orange;

    % Increase counter
    count = count + 1;

    % Specify midesction pattern
    triangles = zeros((2*n)+1, 2);

    % Populating triangles

    % Top sections
    for i = 2:2:2*n
        triangles(i, 1) = ((i-2)/2)*ls + x;
        triangles(i, 2) = h1 + l;
    end

    % Bottom sections
    for j = 1:2:((2*(n)+1))
        triangles(j, 1) = ((j-1)/2)*ls;
        triangles(j, 2) = h1;    
    end

    if triangles(2, 1) == 0

        % Log data to structure
        dataFoldB(count).x = triangles(2:end, 1);
        dataFoldB(count).y = triangles(2:end, 2);
        dataFoldB(count).color = blue;

    else
        % Log data to structure
        dataFoldB(count).x = triangles(:, 1);
        dataFoldB(count).y = triangles(:, 2);
        dataFoldB(count).color = blue;

    end
    
end

if x < 10^-4 || x > ls - 10^-4
    
    % Initialize
    middletube = zeros(2*n, 2);
    
    % Define middle section tube boundary lines
    for i = 0:2:(2*n)

        % Increase count initially
        count = count + 1;

        index = ((i)/2) + 1;
        middletube(i+1, 1) = (index-1)*ls;
        middletube(i+1, 2) = h1;
        middletube(i+2, 1) = (index-1)*ls;
        middletube(i+2, 2) = h1+l;

        % Log data to structure and add to plot. Plotting is sequential
        dataFoldB(count).x = middletube(i+1:i+2, 1);
        dataFoldB(count).y = middletube(i+1:i+2, 2);
        dataFoldB(count).color = blue;

    end  
    
end


% axis off

% Specifying top and bottom tube patterns
%-------------------------------------------

% Bottom tube folds and graphing
bottomtube = zeros(2*(n-1), 2);

% Since we are assuming no offset on the initial triangle from bottom tube,
% we only need 5 side lengths to be sketched (+ 1 additional for connect)
for ii = 0:2:(2*n)-4
    
    % Increase count initially
    count = count + 1;
    
    index = ((ii)/2) + 1;
    bottomtube(ii+1, 1) = index*ls;
    bottomtube(ii+1, 2) = 0;
    bottomtube(ii+2, 1) = index*ls;
    bottomtube(ii+2, 2) = h1;
    
    % Log data to structure and add to plot. Plotting is sequential
    dataFoldB(count).x = bottomtube(ii+1:ii+2, 1);
    dataFoldB(count).y = bottomtube(ii+1:ii+2, 2);
    dataFoldB(count).color = blue;
    
end

% disp(bottomtube)

% Top tube folds and graphing. Provided for x = 0 and x =/= 0 (just in
% case), although it seems like x will always be greater than 0
if x == 0
    toptube = zeros(2*(n-1), 2);
    for jj = 0:2:2*(n-1)-2
        
        % Increase count initially
        count = count + 1;
        
        index = (jj/2) + 1;
        toptube(jj+1, 1) = index*ls;
        toptube(jj+1, 2) = h1 + l;
        toptube(jj+2, 1) = index*ls;
        toptube(jj+2, 2) = h1 + l + h2;
        
        % Log data to structure and add to plot. Plotting is sequential
        dataFoldB(count).x = toptube(jj+1:jj+2, 1);
        dataFoldB(count).y = toptube(jj+1:jj+2, 2);
        dataFoldB(count).color = blue;

    end
    
else
    toptube = zeros(2*(n), 2);
    for jj = 0:2:2*n-2
        
        % Increase count initially
        count = count + 1;
        
        index = (jj/2);
        toptube(jj+1, 1) = index*ls + x;
        toptube(jj+1, 2) = h1 + l;
        toptube(jj+2, 1) = index*ls + x;
        toptube(jj+2, 2) = h1 + l + h2;
        
        % Log data to structure and add to plot. Plotting is sequential
        dataFoldB(count).x = toptube(jj+1:jj+2, 1);
        dataFoldB(count).y = toptube(jj+1:jj+2, 2);
        dataFoldB(count).color = blue;
       
    end
    
end

% Label the plot for clarity
title({
    ('Origami Schematic B for Provided Parameters:')
    ['[r = ' num2str(r) ', n = ' num2str(n) ', h = ' num2str(h) ', alpha = ' num2str(alpha) ']']
    })

daspect([1 1 1])

% m = r*alpha;
m = floor(alpha/((2*pi)/n))*ls + x;
lmax = h1 + l + h2;


end
