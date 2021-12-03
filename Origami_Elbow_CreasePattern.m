% Graph for crease pattern - Origami elbow fitting
% Last edited 6/8/2021 by Lucien Peach

function [dataFoldA, m, lmax] = Origami_Elbow_CreasePattern(lengths, ls, n, h1, h2, r, phi, theta, mirror, split, tuckangle)

% Create "duplicate" value
duplicate = 1;

% Define theta_original
theta_original = theta;

% Check value of theta
if theta > pi/2 && strcmp(split, 'off') ~= 1
    
    duplicate = 2;
    
end

% Counter used for data structure indexing
count = 1;

% Identify colors
orange = [1, 0.41, 0];
% red = [1, 0, 0];
blue = [0, 0, 1];
black = [0, 0, 0];

% Specifying boundary for sheet - in 4 parts
% -------------------------------------------------------------------

% Determine max value of lengths array
lmax = max(lengths);

% Determine if the layer needs to be duplicated or not
if duplicate == 1
    
% Overall boundary
    boundarybottom = [0, 0; n*ls, 0];
    boundaryleft = [0, 0; 0, h1 + 2*lmax + h2];
    boundarytop = [0, h1 + 2*lmax + h2; n*ls, h1 + 2*lmax + h2];
    boundaryright = [n*ls, h1 + 2*lmax + h2; n*ls, 0];

else
    
    % Overall boundary
    boundarybottom = [0, 0; n*ls, 0];
    boundaryleft = [0, 0; 0, h1 + 4*lmax + h2];
    boundarytop = [0, h1 + 4*lmax + h2; n*ls, h1 + 4*lmax + h2];
    boundaryright = [n*ls, h1 + 4*lmax + h2; n*ls, 0];
    
end

% Log data to structure
dataFoldA(count).x = boundarybottom(:, 1);
dataFoldA(count).y = boundarybottom(:, 2);
dataFoldA(count).color = black;

count = count + 1;
dataFoldA(count).x = boundaryleft(:, 1);
dataFoldA(count).y = boundaryleft(:, 2);
dataFoldA(count).color = black;

count = count + 1;
dataFoldA(count).x = boundarytop(:, 1);
dataFoldA(count).y = boundarytop(:, 2);
dataFoldA(count).color = black;

count = count + 1;
dataFoldA(count).x = boundaryright(:, 1);
dataFoldA(count).y = boundaryright(:, 2);
dataFoldA(count).color = black;

% Increase count
count = count + 1;

dataFoldA(count).x = boundaryleft(:, 1);
dataFoldA(count).y = boundaryleft(:, 2);
dataFoldA(count).color = blue;

count = count + 1;

dataFoldA(count).x = boundaryright(:, 1);
dataFoldA(count).y = boundaryright(:, 2);
dataFoldA(count).color = blue;

% Begin by specifying proximal and distal horizontal lines
% -------------------------------------------------------------------

% % Store proximal to array
% proximalx = [0; n*ls];
% proximaly = [h1; h1];
% dataFoldA(count).x = proximalx;
% dataFoldA(count).y = proximaly;
% dataFoldA(count).color = red;
% 
% % Begin plot. Plot proximal line
% plot(dataFoldA(count).x, dataFoldA(count).y, 'color', ...
%     dataFoldA(count).color);
% hold on
% set(gcf, 'color', 'w');
% 
% % Increase counter
% count = count + 1;
% 
% % Store distal to array
% distalx = [0; n*ls];
% distaly = [h1 + 2*lmax; h1 + 2*lmax];
% dataFoldA(count).x = distalx;
% dataFoldA(count).y = distaly;
% dataFoldA(count).color = red;
% 
% % Plot distal line
% plot(dataFoldA(count).x, dataFoldA(count).y, 'color', ...
%     dataFoldA(count).color);
% 
% % Increase count
% count = count + 1;

% Specify line pattern created by lengths
% ------------------------------------------------------------------

% Increase count
count = count + 1;

if duplicate == 1

    % Initialize midsection array
    midsection = zeros(n+1, 2);

    % Store each length value to array. Loop.

    % X coordinates
    for i = 1:n+1
        midsection(i, 1) = (i-1) * ls;    
    end

    % Y coordinates
    for j = 1:n+1
        if j < n+1
            midsection(j, 2) = lengths(j, 1) + h1;
        else
            midsection(j, 2) = midsection(1, 2);
        end  
    end

    % % Log data for final point (extra column for folding)
    % midsection(n+2, 1) = (n+1)*ls;
    % midsection(n+2, 2) = midsection(2, 2);
    
else
    
    % Initialize duplication midsection
    midsection = zeros(2*(n+1), 2);

    % Store each length value to array. Loop.

    % X coordinates
    for i = 1:n+1
        midsection(i, 1) = (i-1) * ls;    
    end
    
    midsection(n+2:2*(n+1), 1) = midsection(1:n+1, 1);

    % Y coordinates
    for j = 1:n+1
        if j < n+1
            midsection(j, 2) = lengths(j, 1) + h1;
        else
            midsection(j, 2) = midsection(1, 2);
        end  
    end
    
    midsection(n+2:2*(n+1), 2) = midsection(1:n+1, 2) + 2*lmax;
end

% Log data to structure
dataFoldA(count).x = midsection(1:n+1, 1);
dataFoldA(count).y = midsection(1:n+1, 2);
dataFoldA(count).color = blue;

% Counter
count = count + 1;

% Log data to structure
dataFoldA(count).x = midsection(n+2:end, 1);
dataFoldA(count).y = midsection(n+2:end, 2);
dataFoldA(count).color = blue;

% Mirror (Optional)
% ------------------------------------------------------------------

if strcmp(mirror, 'on') == 1
    
    if duplicate == 1

        % Counter
        count = count + 1;

        % Mirror lengths region
        mirrormid = zeros(n+1, 2);

        % Store each length value to array. Loop.

        % X coordinates
        for i = 1:n+1
            mirrormid(i, 1) = (i-1) * ls;    
        end

        % Y coordinates
        for j = 1:n+1
            if j < n+1
                mirrormid(j, 2) = -lengths(j, 1) + h1 + 2*lmax;
            else
                mirrormid(j, 2) = mirrormid(1, 2);
            end  
        end
        
        % Log data to structure
        dataFoldA(count).x = mirrormid(:, 1);
        dataFoldA(count).y = mirrormid(:, 2);
        dataFoldA(count).color = blue;
        
    else
        
        % Counter
        count = count + 1;

        % Mirror lengths region
        mirrormid = zeros(2*(n+1), 2);

        % Store each length value to array. Loop.

        % X coordinates
        for i = 1:n+1
            mirrormid(i, 1) = (i-1) * ls;    
        end
        
        mirrormid(n+2:2*(n+1), 1) = mirrormid(1:n+1, 1);

        % Y coordinates
        for j = 1:n+1
            if j < n+1
                mirrormid(j, 2) = -lengths(j, 1) + h1 + 2*lmax;
            else
                mirrormid(j, 2) = mirrormid(1, 2);
            end      
        end
        
        mirrormid(n+2:2*(n+1), 2) = mirrormid(1:n+1, 2) + 2*lmax;
        
        % Log data to structure
        dataFoldA(count).x = mirrormid(1:n+1, 1);
        dataFoldA(count).y = mirrormid(1:n+1, 2);
        dataFoldA(count).color = blue;
        
        % Increase count
        count = count + 1;
        
        dataFoldA(count).x = mirrormid(n+2:2*(n+1), 1);
        dataFoldA(count).y = mirrormid(n+2:2*(n+1), 2);
        dataFoldA(count).color = blue;
        
    end
   
end

% Specify boundary of hidden section
% ------------------------------------------------------------------

% Increase count
count = count + 1;

if duplicate == 1

    % Store hidden to array
    hiddenx = [0; n*ls];
    hiddeny = [h1 + lmax; h1 + lmax];
    dataFoldA(count).x = hiddenx;
    dataFoldA(count).y = hiddeny;
    dataFoldA(count).color = orange;
    
else

    % Store hidden to array
    hiddenx = [0; n*ls];
    hiddeny = [h1 + lmax; h1 + lmax];
    dataFoldA(count).x = hiddenx;
    dataFoldA(count).y = hiddeny;
    dataFoldA(count).color = orange;
    
    % Increase count
    count = count + 1;
    
    % Store hidden2 to array
    hiddenx = [0; n*ls];
    hiddeny = [h1 + 3*lmax; h1 + 3*lmax];
    dataFoldA(count).x = hiddenx;
    dataFoldA(count).y = hiddeny;
    dataFoldA(count).color = orange;
    
end
    
% Specify vertical segments of hidden section
% ------------------------------------------------------------------

if duplicate == 1
    
    % Initialize hidden vertical line array
    hiddenmid = zeros(2*(n-1), 2);

    % Loop through to store [x,y] pairs for each line segment
    for k = 0:2:(2*n)-4

        % Increase count initially 
        count = count + 1;

        % Indexing 
        index = (k/2) + 2;

        % Store x coordinates for each segment
        hiddenmid(k+1, 1) = (index-1) * ls;
        hiddenmid(k+2, 1) = (index-1) * ls;

        % Store y coordinates for each segment
        hiddenmid(k+1, 2) = lengths(index, 1) + h1;
        hiddenmid(k+2, 2) = h1 + lmax;

        % Store to data array
        dataFoldA(count).x = hiddenmid(k+1:k+2, 1);
        dataFoldA(count).y = hiddenmid(k+1:k+2, 2);
        dataFoldA(count).color = orange;

    end
    
else
    
    % Initialize hidden vertical line array
    hiddenmid = zeros(2*(2*(n-1)), 2);

    % Loop through to store [x,y] pairs for each line segment
    for k = 0:2:(2*n)-4

        % Increase count initially 
        count = count + 1;

        % Indexing 
        index = (k/2) + 2;

        % Store x coordinates for each segment
        hiddenmid(k+1, 1) = (index-1) * ls;
        hiddenmid(k+2, 1) = (index-1) * ls;

        % Store y coordinates for each segment
        hiddenmid(k+1, 2) = lengths(index, 1) + h1;
        hiddenmid(k+2, 2) = h1 + lmax;

        % Store to data array
        dataFoldA(count).x = hiddenmid(k+1:k+2, 1);
        dataFoldA(count).y = hiddenmid(k+1:k+2, 2);
        dataFoldA(count).color = orange;

    end
    
    hiddenmid(2*(n-1)+1:end, 1) = hiddenmid(1:2*(n-1), 1);
    hiddenmid(2*(n-1)+1:end, 2) = hiddenmid(1:2*(n-1), 2) + 2*lmax;
    
    for i = 2*(n-1)+1:2:4*(n-1)-1
        
        % counter
        count = count + 1;
        
        % Store to data array
        dataFoldA(count).x = hiddenmid(i:i+1, 1);
        dataFoldA(count).y = hiddenmid(i:i+1, 2);
        dataFoldA(count).color = orange;
        
    end
    
end

% Bottom tube folds
% --------------------------------------------------------------------

if duplicate == 1
    
    % Bottom tube folds and graphing
    bottomtube = zeros(2*(n-1), 2);

    % Ignore side folds as these are graphed by the boundary section
    for ii = 0:2:(2*n)-4

        % Increase count initially
        count = count + 1;

        % Indexing 
        index = ((ii)/2) + 2;

        % Populate array
        bottomtube(ii+1, 1) = (index-1)*ls;
        bottomtube(ii+1, 2) = 0;
        bottomtube(ii+2, 1) = (index-1)*ls;
        bottomtube(ii+2, 2) = lengths(index, 1) + h1;

        % Log data to structure and add to plot. Plotting is sequential
        dataFoldA(count).x = bottomtube(ii+1:ii+2, 1);
        dataFoldA(count).y = bottomtube(ii+1:ii+2, 2);
        dataFoldA(count).color = blue;

    end

else
    
    % Bottom tube folds and graphing
    bottomtube = zeros(2*(2*(n-1)), 2);

    % Ignore side folds as these are graphed by the boundary section
    for ii = 0:2:(2*n)-4

        % Increase count initially
        count = count + 1;

        % Indexing 
        index = ((ii)/2) + 2;

        % Populate array
        bottomtube(ii+1, 1) = (index-1)*ls;
        bottomtube(ii+1, 2) = 0;
        bottomtube(ii+2, 1) = (index-1)*ls;
        bottomtube(ii+2, 2) = lengths(index, 1) + h1;

        % Log data to structure and add to plot. Plotting is sequential
        dataFoldA(count).x = bottomtube(ii+1:ii+2, 1);
        dataFoldA(count).y = bottomtube(ii+1:ii+2, 2);
        dataFoldA(count).color = blue;

    end
    
    % For duplication
    bottomtube(2*(n-1)+1:end, 1) = bottomtube(1:2*(n-1), 1);
    bottomtube(2*(n-1)+1:2:end, 2) = bottomtube(1:2:2*(n-1), 2) + 2*lmax + h1;
    bottomtube(2*(n-1)+2:2:end, 2) = bottomtube(2:2:2*(n-1), 2) + 2*lmax;
    
    % Ignore side folds as these are graphed by the boundary section
    for i = 2*(n-1)+1:2:2*(2*(n-1))-1
        
        % counter
        count = count + 1;

        % Log data to structure and add to plot. Plotting is sequential
        dataFoldA(count).x = bottomtube(i:i+1, 1);
        dataFoldA(count).y = bottomtube(i:i+1, 2);
        dataFoldA(count).color = blue;

    end
    
end

% Top tube folds
% --------------------------------------------------------------------

if duplicate == 1
    
    % Top tube folds and graphing
    toptube = zeros(2*(n-1), 2);

    % Ignore side folds as these are graphed by the boundary section
    for jj = 0:2:(2*n)-4

        % Increase count initially
        count = count + 1;

        % Indexing 
        index = (jj/2) + 2;

        % Populate array
        toptube(jj+1, 1) = (index-1)*ls;
        toptube(jj+1, 2) = lmax + h1;
        toptube(jj+2, 1) = (index-1)*ls;
        toptube(jj+2, 2) = 2*lmax + h1 + h2;

        % Log data to structure and add to plot. Plotting is sequential
        dataFoldA(count).x = toptube(jj+1:jj+2, 1);
        dataFoldA(count).y = toptube(jj+1:jj+2, 2);
        dataFoldA(count).color = blue;

    end


else
    
    % Top tube folds and graphing
    toptube = zeros(2*(2*(n-1)), 2);

    % Ignore side folds as these are graphed by the boundary section
    for jj = 0:2:(2*n)-4

        % Increase count initially
        count = count + 1;

        % Indexing 
        index = (jj/2) + 2;

        % Populate array
        toptube(jj+1, 1) = (index-1)*ls;
        toptube(jj+1, 2) = lmax + h1;
        toptube(jj+2, 1) = (index-1)*ls;
        toptube(jj+2, 2) = 2*lmax + h1;

        % Log data to structure and add to plot. Plotting is sequential
        dataFoldA(count).x = toptube(jj+1:jj+2, 1);
        dataFoldA(count).y = toptube(jj+1:jj+2, 2);
        dataFoldA(count).color = blue;
        
    end
    
    % Add duplication data
    toptube(2*(n-1)+1:end, 1) = toptube(1:2*(n-1), 1);
    toptube(2*(n-1)+1:2:end, 2) = toptube(1:2:2*(n-1), 2) + 2*lmax;
    toptube(2*(n-1)+2:2:end, 2) = toptube(2:2:2*(n-1), 2) + 2*lmax + h2;
    
    % Loop to plot
    for i = 2*(n-1)+1:2:2*(2*(n-1))-1
               
        % counter
        count = count + 1;

        % Log data to structure and add to plot. Plotting is sequential
        dataFoldA(count).x = toptube(i:i+1, 1);
        dataFoldA(count).y = toptube(i:i+1, 2);
        dataFoldA(count).color = blue;

    end
    
end

% Tuck Angle Additions
% -----------------------------------------------------------------------

% First, find the tucklengths vector that will contain the tucklength
% values, calculated using the value of max(lengths) as well as the
% individual length segments.
tucklengths = zeros(n+1, 1);

for i = 1:n+1
    
    tucklengths(i, 1) = max(lengths) - lengths(i, 1);
    
end

% Using this value as well as the values contained within tuckangle, we can
% determine the "hidden region" of the schematic
tuckoffsets = zeros(n, 1);

for i = 1:n
    
    tuckoffsets(i, 1) = tan(tuckangle(i, 1))*tucklengths(i, 1);
    
end

% We can now use this array to print the overlap region for the print
% schematic for the elbow joint
tuckarray = zeros(3*n, 2);

% Also create "small" array for the instance of no mirroring
tuckarray_small = zeros(2*n, 2);

% There are 4 cases to address here.
% Note that the execution of Case 1 must precede the execution of Case 4,
% and that the execution of Case 2 must precede the execution of Case 3.
% This is why the duplicate value is not included in the conditional
% statements for Case 1 and Case 2.

% Case 1: mirror = 'on'
if strcmp(mirror, 'on') == 1
    for i = 1:n

        % Array index initialization
        index = (i-1)*3 + 1;

        % Increase count
        count = count + 1;

        % Storing plotting points to array
        tuckarray(index, 1) = (i-1)*ls;
        tuckarray(index, 2) = lengths(i, 1) + h1;

        tuckarray(index+1, 1) = ((i-1)*ls) + tuckoffsets(i);
        tuckarray(index+1, 2) = max(lengths) + h1;

        tuckarray(index+2, 1) = (i-1)*ls;
        tuckarray(index+2, 2) = 2*max(lengths) - lengths(i, 1) + h1;

        % Storing to DataStruct and plotting
        dataFoldA(count).x = tuckarray(index:index+2, 1);
        dataFoldA(count).y = tuckarray(index:index+2, 2);
        dataFoldA(count).color = orange;

    end

% Case 2: duplicate = 1, mirror = 'off'
elseif strcmp(mirror, 'off') == 1
    for i = 1:n

        % Array index initialization
        index = (i-1)*2 + 1;

        % Increase count
        count = count + 1;

        % Storing plotting points to array
        tuckarray_small(index, 1) = (i-1)*ls;
        tuckarray_small(index, 2) = lengths(i, 1) + h1;

        tuckarray_small(index+1, 1) = ((i-1)*ls) + tuckoffsets(i);
        tuckarray_small(index+1, 2) = max(lengths) + h1;

        % Storing to DataStruct and plotting
        dataFoldA(count).x = tuckarray_small(index:index+1, 1);
        dataFoldA(count).y = tuckarray_small(index:index+1, 2);
        dataFoldA(count).color = orange;

    end
    
end
        
% Case 3: duplicate = 2, mirror = 'off'        
if duplicate == 2 && strcmp(mirror, 'off') == 1
% Create the array for the secondary level, if needed, by adding
% 2*max(lengths) to all values of the primary level
   
    % Initialize array for small array (no mirroring)
    newtuckarray_small = zeros(size(tuckarray_small, 1), 2);
    
    % Populate new array 
    for i = 1:size(tuckarray_small, 1)
    
        newtuckarray_small(i, 1) = tuckarray_small(i, 1);
        newtuckarray_small(i, 2) = tuckarray_small(i, 2) + 2*max(lengths);
    
    end
    
    % Store array data to structure and plot
    for i = 1:n
        
        count = count + 1;
        
        index = (i-1)*2 + 1;
        
        dataFoldA(count).x = newtuckarray_small(index:index+1, 1);
        dataFoldA(count).y = newtuckarray_small(index:index+1, 2);
        dataFoldA(count).color = orange;
        
    end
    
% Case 4: duplicate = 2, mirror = 'on'
elseif duplicate == 2 && strcmp(mirror, 'on') == 1
% Create the array for the secondary level, if needed, by adding
% 2*max(lengths) to all values of the primary level
   
    % Initialize array
    newtuckarray = zeros(size(tuckarray, 1), 2);
    
    % Populate new array 
    for i = 1:size(tuckarray, 1)
    
        newtuckarray(i, 1) = tuckarray(i, 1);
        newtuckarray(i, 2) = tuckarray(i, 2) + 2*max(lengths);
    
    end
    
    % Store array data to structure and plot
    for i = 1:n
        
        count = count + 1;
        
        index = (i-1)*3 + 1;
        
        dataFoldA(count).x = newtuckarray(index:index+2, 1);
        dataFoldA(count).y = newtuckarray(index:index+2, 2);
        dataFoldA(count).color = orange;
        
    end
     
end

% Plotting Options
% ----------------------------------------------------------------------

% Label the plot for clarity
% title({
%     ('Origami Schematic A for Provided Parameters:')
%     ['[r = ' num2str(r) ', n = ' num2str(n) ', phi = ' num2str(phi) ', theta = ' num2str(theta_original) ']']
%     })
% 
% daspect([1 1 1])

m = 0;

if duplicate == 1
    lmax = h1 + 2*lmax + h2;
else
    lmax = h1 + 4*lmax + h2;
end

end
