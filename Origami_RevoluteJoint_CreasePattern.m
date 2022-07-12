function [dataFoldD, m, lmax] = Origami_RevoluteJoint_CreasePattern(lengths, ls, n, h1, h2, nz)
% ORIGAMI_REVOLUTEJOINT_CREASEPATTERN - Generates a crease pattern for the
% origami REBO joint.

% Inputs:
%   lengths     - vector of length measurements for use in plotting
%                 revolute joint section of schematic. 
%   ls          - side length of folded origami linkage. 
%   n           - number of sides of folded origami linkage. 
%   h1          - length from bottom of schematic to base of revolute
%                 section. ie height of lower tube region. 
%   h2          - length from upper boundary of revolute section to top of
%                 full schematic. ie height of upper tube region. 
%   nz          - number of recursive sink layers.

% Outputs:
%   dataFoldD   - data structure which contains pertinent information for
%                 DXF generation of crease schematic. 
%   m           - horizontal offset value for use within DataFoldAppend.m.
%   lmax        - total height of crease schematic. Will be used in crease
%                 "stacking" and duplication. 

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last edited 6/11/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Counter used for data structure indexing
count = 1;

% Identify colors
orange = [1, 0.41, 0];
blue = [0, 0, 1];
% orange = [0, 0, 1];
% blue = [1, 0, 0];
black = [0, 0, 0];
% yellow = [1, 1, 0];

% Begin by specifying sheet boundary
% -------------------------------------------------------------------

% Determine max value of lengths array
lmax = max(lengths);

% Store boundary coordinates to array
boundarybottom = [0, 0; n*ls, 0];
boundaryleft = [0, 0; 0, h1 + h2 + 2*lmax];
boundarytop = [0, h1 + h2 + 2*lmax; n*ls, h1 + h2 + 2*lmax];
boundaryright = [n*ls, h1 + h2 + 2*lmax; n*ls, 0];

dataFoldD(count).x = boundarybottom(:, 1);
dataFoldD(count).y = boundarybottom(:, 2);
dataFoldD(count).color = black;

count = count + 1;
dataFoldD(count).x = boundaryleft(:, 1);
dataFoldD(count).y = boundaryleft(:, 2);
dataFoldD(count).color = black;

count = count + 1;
dataFoldD(count).x = boundarytop(:, 1);
dataFoldD(count).y = boundarytop(:, 2);
dataFoldD(count).color = black;

count = count + 1;
dataFoldD(count).x = boundaryright(:, 1);
dataFoldD(count).y = boundaryright(:, 2);
dataFoldD(count).color = black;

% Increase counter
count = count + 1;

dataFoldD(count).x = boundaryleft(:, 1);
dataFoldD(count).y = boundaryleft(:, 2);
dataFoldD(count).color = blue;

count = count + 1;
dataFoldD(count).x = boundaryright(:, 1);
dataFoldD(count).y = boundaryright(:, 2);
dataFoldD(count).color = blue;

% Increase counter
count = count + 1;

% Specify horizontal red lines (proximal / distal)
% ------------------------------------------------------------------

% Proximal line
proximal = [0, h1; n*ls, h1];
dataFoldD(count).x = proximal(:, 1);
dataFoldD(count).y = proximal(:, 2);
dataFoldD(count).color = blue;

% Increase counter
count = count + 1;

% Distal line
distal = [0, h1 + (2*lmax); n*ls, h1 + (2*lmax)];
dataFoldD(count).x = distal(:, 1);
dataFoldD(count).y = distal(:, 2);
dataFoldD(count).color = blue;

% Bottom tube folds
% --------------------------------------------------------------------

% Bottom tube folds and graphing
bottomtube = zeros(2*n, 2);

% Ignore side folds as these are graphed by the boundary section
for ii = 1:2:2*(n-1)
    
    % Increase count initially
    count = count + 1;
    
    % Indexing 
    index = ((ii-1)/2) + 1;
    
    % Populate array
    bottomtube(ii, 1) = index*ls;
    bottomtube(ii, 2) = 0;
    bottomtube(ii+1, 1) = index*ls;
    bottomtube(ii+1, 2) = h1;
    
    % Log data to structure and add to plot. Plotting is sequential
    dataFoldD(count).x = bottomtube(ii:ii+1, 1);
    dataFoldD(count).y = bottomtube(ii:ii+1, 2);
    dataFoldD(count).color = blue;

end

% Top tube folds
% --------------------------------------------------------------------

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
    toptube(jj, 2) = (2*lmax) + h1;
    toptube(jj+1, 1) = index*ls;
    toptube(jj+1, 2) = (2*lmax) + h1 + h2;
    
    % Log data to structure and add to plot. Plotting is sequential
    dataFoldD(count).x = toptube(jj:jj+1, 1);
    dataFoldD(count).y = toptube(jj:jj+1, 2);
    dataFoldD(count).color = blue;

end

% Cross hatch pattern
% --------------------------------------------------------------------

% The number of cross hatches will always be n-2. We can use this
% information to initialize the cross hatch array. Each cross hatch
% requires 4 coordinate points
cross_hatch = zeros(4*(n-1), 2);

% Populating cross_hatch will be a delicate task. We will use a loop inside
% a loop to account for the 2x repetitive nature of this pattern

for i = 1:2
    for j = 0:4:(2*(n-1) - 4) % Accounts for each 1/2 segment of total
        if i == 1
            
            % Indexing for initial half of cross hatches
            index = (j/4) + 1;
            
            % Loop again for increasing / decreasing hatches
            for k = 1:2
                
                % Increase count initially
                count = count + 1;
                if k == 1
                    
                    % Create cross hatch pattern for first batch (X coord)
                    cross_hatch((j+1), 1) = (index-1)*ls;
                    cross_hatch((j+2), 1) = index*ls;
                    
                    % Create cross hatch pattern for first batch (Y coord)
                    cross_hatch((j+1), 2) = h1 + 2*lmax;
                    cross_hatch((j+2), 2) = h1;
                    
                    % Express in structure format and plot iteratively
                    dataFoldD(count).x = cross_hatch(j+1:j+2, 1);
                    dataFoldD(count).y = cross_hatch(j+1:j+2, 2);
                    dataFoldD(count).color = orange;
                    
                else
                    
                    % Create cross hatch pattern for second batch (X coord)
                    cross_hatch((j+3), 1) = index*ls;
                    cross_hatch((j+4), 1) = (index-1)*ls;
                    
                    % Create cross hatch pattern for second batch (Y coord)
                    cross_hatch((j+3), 2) = h1 + 2*lmax;
                    cross_hatch((j+4), 2) = h1;
                    
                    % Express in structure format and plot iteratively
                    dataFoldD(count).x = cross_hatch(j+3:j+4, 1);
                    dataFoldD(count).y = cross_hatch(j+3:j+4, 2);
                    dataFoldD(count).color = orange;
                    
                end
            end
            
            
        else
            
            % Shifts to second regime
            offset = 2*(n-2);
            
            % Indexing for second half of cross hatches
            index2 = (j/4) + (n/2) + 1;
            
            % Loop again for increasing / decreasing hatches
            for kk = 1:2
                
                % Increase count initially
                count = count + 1;
                if kk == 1
                    
                    % Create cross hatch pattern for first batch (X coord)
                    cross_hatch((j+1+offset), 1) = (index2-1)*ls;
                    cross_hatch((j+2+offset), 1) = index2*ls;
                    
                    % Create cross hatch pattern for first batch (Y coord)
                    cross_hatch((j+1+offset), 2) = h1 + 2*lmax;
                    cross_hatch((j+2+offset), 2) = h1;
                    
                    % Express in structure format and plot iteratively
                    dataFoldD(count).x = cross_hatch(j+1+offset:j+2+offset, 1);
                    dataFoldD(count).y = cross_hatch(j+1+offset:j+2+offset, 2);
                    dataFoldD(count).color = orange;
                    
                else
                    
                    % Create cross hatch pattern for second batch (X coord)
                    cross_hatch((j+3+offset), 1) = index2*ls;
                    cross_hatch((j+4+offset), 1) = (index2-1)*ls;
                    
                    % Create cross hatch pattern for second batch (Y coord)
                    cross_hatch((j+3+offset), 2) = h1 + 2*lmax;
                    cross_hatch((j+4+offset), 2) = h1;
                    
                    % Express in structure format and plot iteratively
                    dataFoldD(count).x = cross_hatch(j+3+offset:j+4+offset, 1);
                    dataFoldD(count).y = cross_hatch(j+3+offset:j+4+offset, 2);
                    dataFoldD(count).color = orange;
                    
                end
            end            
        end 
    end
end

% % Append final section for glue-over segment
% count = count + 1;
% 
% % Assign to cross_hatch with a horizontal shift
% cross_hatch(end-3, 1) = cross_hatch(1, 1) + n*ls;
% cross_hatch(end-3, 2) = cross_hatch(1, 2);
% cross_hatch(end-2, 1) = cross_hatch(2, 1) + n*ls;
% cross_hatch(end-2, 2) = cross_hatch(2, 2);
% 
% % Storing to array
% dataFoldD(count).x = cross_hatch(end-3:end-2, 1);
% dataFoldD(count).y = cross_hatch(end-3:end-2, 2);
% dataFoldD(count).color = orange;
% 
% % Plotting
% plot(dataFoldD(count).x, dataFoldD(count).y, 'color', ...
%     dataFoldD(count).color);
% 
% % Repeat for second cross hatch
% count = count + 1;
% 
% % Assign to cross_hatch with a horizontal shift
% cross_hatch(end-1, 1) = cross_hatch(3, 1) + n*ls;
% cross_hatch(end-1, 2) = cross_hatch(3, 2);
% cross_hatch(end, 1) = cross_hatch(4, 1) + n*ls;
% cross_hatch(end, 2) = cross_hatch(4, 2);
% 
% % Storing to array
% dataFoldD(count).x = cross_hatch(end-1:end, 1);
% dataFoldD(count).y = cross_hatch(end-1:end, 2);
% dataFoldD(count).color = orange;
% 
% % Plotting
% plot(dataFoldD(count).x, dataFoldD(count).y, 'color', ...
%     dataFoldD(count).color);


% Create "weave" pattern (pt.1) - horizontal segments (orange)
% ------------------------------------------------------------------

% Initialization of orange horizontal array
o_val = (n/2)-2;
weave2 = zeros((o_val*4) + 4, 2);

% Increase count
count = count + 1;

% End segment orange horizontal, first set
weave2((o_val*4)+1, 1) = ((n/2)-1)*ls;
weave2((o_val*4)+2, 1) = (n/2)*ls;
weave2((o_val*4)+1, 2) = h1 + lengths(1);
weave2((o_val*4)+2, 2) = h1 + lengths(1);

dataFoldD(count).x = weave2((o_val*4)+1:(o_val*4)+2, 1);
dataFoldD(count).y = weave2((o_val*4)+1:(o_val*4)+2, 2);
dataFoldD(count).color = orange;

% Increase count
count = count + 1;

% End segment orange horizontal, second set
weave2((o_val*4)+3, 1) = (n-1)*ls;
weave2((o_val*4)+4, 1) = n*ls;
weave2((o_val*4)+3, 2) = h1 + lengths(1);
weave2((o_val*4)+4, 2) = h1 + lengths(1);

dataFoldD(count).x = weave2((o_val*4)+3:(o_val*4)+4, 1);
dataFoldD(count).y = weave2((o_val*4)+3:(o_val*4)+4, 2);
dataFoldD(count).color = orange;

% Looping for horizontal segments within diamond regions
% Segment 1
for ii = 1:2:(o_val*2)-1
    if n > 4 % Terms n=4 and less do not have this component

        % Increase count initially
        count = count + 1;

        % Indexing
        index = (ii-1)/2 + 1;

        % Populating weave2
        weave2(ii, 1) = (index-1)*ls + ls/2;
        weave2(ii, 2) = h1 + lengths(1);
        weave2(ii+1, 1) = (index)*ls + ls/2;
        weave2(ii+1, 2) = h1 + lengths(1);

        % Store values to structure
        dataFoldD(count).x = weave2(ii:ii+1, 1);
        dataFoldD(count).y = weave2(ii:ii+1, 2);
        dataFoldD(count).color = orange;
    
    end
    
end

% Segment 2
for ii = (o_val*2):2:(o_val*4)-2
    if n > 4 % Terms n=4 and less do not have this component

        % Increase count initially
        count = count + 1;

        % Indexing
        index = (ii-(o_val*2))/2 + 1 + n/2;

        % Populating weave2
        weave2(ii+1, 1) = (index-1)*ls + ls/2;
        weave2(ii+1, 2) = h1 + lengths(1);
        weave2(ii+2, 1) = (index)*ls + ls/2;
        weave2(ii+2, 2) = h1 + lengths(1);

        % Store values to structure
        dataFoldD(count).x = weave2(ii+1:ii+2, 1);
        dataFoldD(count).y = weave2(ii+1:ii+2, 2);
        dataFoldD(count).color = orange;
        
    end
end


% Create "weave" pattern (pt.2) - horizontal segments (blue)
% ------------------------------------------------------------------

% Initialization
quant = (1 + 2*((n/2) - 2)) + 4;
weave1 = zeros(2*quant, 2);

% Final 8 points constitute the stagnant lines (not within loop)

% Increase count
count = count + 1;

% Begin with first horizontal line for first segment (since breaks 2x pattern)
weave1(end-7, 1) = 0;
weave1(end-6, 1) = ls/2;
weave1(end-7, 2) = h1 + lengths(1);
weave1(end-6, 2) = h1 + lengths(1);

dataFoldD(count).x = weave1(end-7:end-6, 1);
dataFoldD(count).y = weave1(end-7:end-6, 2);
dataFoldD(count).color = blue;

% Increase count
count = count + 1;

% First line horizontal for second segment (since breaks 2x pattern)
weave1(end-3, 1) = (n/2)*ls;
weave1(end-2, 1) = (n/2)*ls + ls/2;
weave1(end-3, 2) = h1 + lengths(1);
weave1(end-2, 2) = h1 + lengths(1);

dataFoldD(count).x = weave1(end-3:end-2, 1);
dataFoldD(count).y = weave1(end-3:end-2, 2);
dataFoldD(count).color = blue;

% Increase count
count = count + 1;

% Final horizontal for segment 1
weave1(end-5, 1) = ((n/2)-1)*ls-(ls/2);
weave1(end-4, 1) = ((n/2)-1)*ls;
weave1(end-5, 2) = h1 + lengths(1);
weave1(end-4, 2) = h1 + lengths(1);

dataFoldD(count).x = weave1(end-5:end-4, 1);
dataFoldD(count).y = weave1(end-5:end-4, 2);
dataFoldD(count).color = blue;

% Increase count
count = count + 1;

% Final horizontal for second segment
weave1(end-1, 1) = (n-1)*ls-(ls/2);
weave1(end, 1) = (n-1)*ls;
weave1(end-1, 2) = h1 + lengths(1);
weave1(end, 2) = h1 + lengths(1);

dataFoldD(count).x = weave1((2*quant)-1:2*quant, 1);
dataFoldD(count).y = weave1((2*quant)-1:2*quant, 2);
dataFoldD(count).color = blue;

% Create "weave" pattern (pt.3) - ridge
% ------------------------------------------------------------------

% In a similar process as before, looping inside looping etc. to populate
% the weave1 array
if n > 4
    for i = 1:2
        for j = 1:(quant-4) % Accounts for each 1/2 segment of total
            if i == 1 % Populating arrays

                % Indexing for initial half of ridges
                if rem(j, 2) == 0 % Even indexes (lower)
                    index = j/2;
                    lengthsindex = (j/2) + 1;
                    weave1(j, 1) = index*ls;
                    weave1(j, 2) = h1 + lengths(lengthsindex);
                else % Odd indexes (upper parts of ridge)
                    index = (j-1)/2; 
                    weave1(j, 1) = index*ls + ls/2;
                    weave1(j, 2) = h1 + max(lengths);
                end


            else

                offset = (quant-4); % For horizontal shift

                % Indexing for second half of ridges
                if rem(j, 2) == 0 % Even indexes (lower)
                    index = (j/2) + (n/2);
                    lengthsindex = (j/2) + (n/2) + 1;
                    weave1(offset+j, 1) = index*ls;
                    weave1(offset+j, 2) = h1 + lengths(lengthsindex);
                else % Odd indexes (upper parts of ridge)
                    index = ((j-1)/2) + (n/2);
                    weave1(offset+j, 1) = index*ls + ls/2;
                    weave1(offset+j, 2) = h1 + max(lengths);
                end

            end
        end

        % Increase count initially
        count = count + 1;

        if i == 1 % Storing to data structure and plotting

            dataFoldD(count).x = weave1(1:(quant-4), 1);
            dataFoldD(count).y = weave1(1:(quant-4), 2);
            dataFoldD(count).color = blue;

        else
            offset2 = (quant-3); % For horizontal shift

            dataFoldD(count).x = weave1(offset2:((2*quant)-8), 1);
            dataFoldD(count).y = weave1(offset2:((2*quant)-8), 2);
            dataFoldD(count).color = blue;

        end
    end
end

% % Add overlap section
% 
% % Initialize 2x2 array
% overlap_ridge = zeros(2, 2);
% 
% % Increase counter
% count = count + 1;
% 
% % Populate array
% overlap_ridge(1, 1) = weave1(1, 1) + n*ls;
% overlap_ridge(1, 2) = weave1(1, 2);
% overlap_ridge(2, 1) = weave1(2, 1) + n*ls;
% overlap_ridge(2, 2) = weave1(2, 2);
% 
% % Store to Data Structure
% dataFoldD(count).x = overlap_ridge(1:2, 1);
% dataFoldD(count).y = overlap_ridge(1:2, 2);
% dataFoldD(count).color = blue;
% 
% % Plot
% plot(dataFoldD(count).x, dataFoldD(count).y, 'color', ...
%     dataFoldD(count).color);

% Create "weave" pattern (pt.4) - vertical lines (blue)
% ------------------------------------------------------------------

% Some of these lines will be present regardless while others will depend
% on the value of n. We can start with the 4 that are unchanged.
vert_blue1 = zeros(8, 2);

% Increase count
count = count + 1;

% Populate array for each segment
vert_blue1(1, 1) = ((n/2)-1)*ls;
vert_blue1(1, 2) = h1;
vert_blue1(2, 1) = ((n/2)-1)*ls;
vert_blue1(2, 2) = h1 + 2*max(lengths);

dataFoldD(count).x = vert_blue1(1:2, 1);
dataFoldD(count).y = vert_blue1(1:2, 2);
dataFoldD(count).color = blue;

% Increase count and proceed to second set (2.1 segment)
count = count + 1;

% Populate array for segment
vert_blue1(3, 1) = (n/2)*ls;
vert_blue1(3, 2) = h1;
vert_blue1(4, 1) = (n/2)*ls;
vert_blue1(4, 2) = h1 + 2*max(lengths);

dataFoldD(count).x = vert_blue1(3:4, 1);
dataFoldD(count).y = vert_blue1(3:4, 2);
dataFoldD(count).color = blue;

% Increase count and proceed to final full segment (2.2 segment)
count = count + 1;

% Populate array for segment
vert_blue1(5, 1) = (n-1)*ls;
vert_blue1(5, 2) = h1;
vert_blue1(6, 1) = (n-1)*ls;
vert_blue1(6, 2) = h1 + 2*max(lengths);

dataFoldD(count).x = vert_blue1(5:6, 1);
dataFoldD(count).y = vert_blue1(5:6, 2);
dataFoldD(count).color = blue;

% Retrace to initial segment for use within overlap patterns
count = count + 1;

% Populate array for segment
vert_blue1(7, 1) = 0;
vert_blue1(7, 2) = h1;
vert_blue1(8, 1) = 0;
vert_blue1(8, 2) = h1 + 2*max(lengths);

dataFoldD(count).x = vert_blue1(7:8, 1);
dataFoldD(count).y = vert_blue1(7:8, 2);
dataFoldD(count).color = blue;

% Now, we need to move on the recurring half segments that are contained
% within each of the cross hatches
if n > 4
    % Initialization
    vert_blue2 = zeros(4*(n-4), 2);
    
    for i = 1:2 % Identifies first and second sections, as before
        if i == 1 % First section
            for j = 0:4:(2*(n-4))-4
                
                % Indexing for horizontal coordinate position
                index = (j/4)+1;
                
                % Populate array with top values
                vert_blue2(j+1, 1) = index*ls;
                vert_blue2(j+1, 2) = max(lengths) + h1;
                vert_blue2(j+2, 1) = index*ls;
                vert_blue2(j+2, 2) = 2*max(lengths) + h1;
                
                % Increase count and plot top values
                count = count + 1;
                dataFoldD(count).x = vert_blue2(j+1:j+2, 1);
                dataFoldD(count).y = vert_blue2(j+1:j+2, 2);
                dataFoldD(count).color = blue;
                
                % Lengths index for varying values
                lengthsindex = (j/4) + 2;
                
                % Populate array with bottom values
                vert_blue2(j+3, 1) = index*ls;
                vert_blue2(j+3, 2) = h1;
                vert_blue2(j+4, 1) = index*ls;
                vert_blue2(j+4, 2) = lengths(lengthsindex) + h1;
                
                % Increase count and plot top values
                count = count + 1;
                dataFoldD(count).x = vert_blue2(j+3:j+4, 1);
                dataFoldD(count).y = vert_blue2(j+3:j+4, 2);
                dataFoldD(count).color = blue;

            end
        else
            for j = 2*(n-4):4:(4*(n-4))-4
                                
                % Indexing for horizontal coordinate position
                index = (n/2) + ((j-2*(n-4))/4) + 1;
                
                % Populate array with top values
                vert_blue2(j+1, 1) = index*ls;
                vert_blue2(j+1, 2) = max(lengths) + h1;
                vert_blue2(j+2, 1) = index*ls;
                vert_blue2(j+2, 2) = 2*max(lengths) + h1;
                
                % Increase count and plot top values
                count = count + 1;
                dataFoldD(count).x = vert_blue2(j+1:j+2, 1);
                dataFoldD(count).y = vert_blue2(j+1:j+2, 2);
                dataFoldD(count).color = blue;
                
                % Lengths index for varying values
                lengthsindex = (n/2) + ((j-2*(n-4))/4) + 2;
                
                % Populate array with bottom values
                vert_blue2(j+3, 1) = index*ls;
                vert_blue2(j+3, 2) = h1;
                vert_blue2(j+4, 1) = index*ls;
                vert_blue2(j+4, 2) = lengths(lengthsindex) + h1;
                
                % Increase count and plot top values
                count = count + 1;
                dataFoldD(count).x = vert_blue2(j+3:j+4, 1);
                dataFoldD(count).y = vert_blue2(j+3:j+4, 2);
                dataFoldD(count).color = blue;

            end
        end
    end
end

% Create "weave" pattern (pt.5) - vertical lines (orange)
% ------------------------------------------------------------------

% Similar to the second-half process from before but much shorter
if n > 4
    % Initialization
    vert_orange = zeros(2*(n-4), 2);
    
    for i = 1:2 % Identifies first or second section
        if i == 1
            for j = 0:2:(n-4)-2
                
                % Indexing for horizontal offset
                index = (j/2) + 1;
                
                % Lengths index
                lengthsindex = (j/2) + 2;
                
                % Populate array
                vert_orange(j+1, 1) = index*ls;
                vert_orange(j+1, 2) = h1 + lengths(lengthsindex);
                vert_orange(j+2, 1) = index*ls;
                vert_orange(j+2, 2) = h1 + max(lengths);
                
                % Increase count, add to structure, and plot
                count = count + 1;
                dataFoldD(count).x = vert_orange(j+1:j+2, 1);
                dataFoldD(count).y = vert_orange(j+1:j+2, 2);
                dataFoldD(count).color = orange;

            end
        else
            for j = (n-4):2:2*(n-4)-2
                
                % Indexing for horizontal offset
                index = (n/2) + ((j-(n-4))/2) + 1;
                
                % Lengths index
                lengthsindex = (n/2) + ((j-(n-4))/2) + 2;
                
                % Populate array
                vert_orange(j+1, 1) = index*ls;
                vert_orange(j+1, 2) = h1 + lengths(lengthsindex);
                vert_orange(j+2, 1) = index*ls;
                vert_orange(j+2, 2) = h1 + max(lengths);
                
                % Increase count, add to structure, and plot
                count = count + 1;
                dataFoldD(count).x = vert_orange(j+1:j+2, 1);
                dataFoldD(count).y = vert_orange(j+1:j+2, 2);
                dataFoldD(count).color = orange;

            end
        end
    end
end

% Create "weave" pattern (pt.7) - top shell (blue)
% ------------------------------------------------------------------

% Again, similar to the second-half process from before
if n > 4
    % Initialization
    blue_shell = zeros(2*(n-4), 2);
    
    for i = 1:2 % Identifies first or second section
        if i == 1
            for j = 0:2:(n-4)-2
                
                % Indexing for horizontal offset
                index = (j/2) + 1;
                
                % Lengths index
                lengthsindex = (j/2) + 2;
                
                % Populate array
                blue_shell(j+1, 1) = index*ls;
                blue_shell(j+1, 2) = h1 + 2*max(lengths) - lengths(lengthsindex);
                blue_shell(j+2, 1) = index*ls + ls/2;
                blue_shell(j+2, 2) = h1 + max(lengths);
                
                % Increase count, add to structure, and plot
                count = count + 1;
                dataFoldD(count).x = blue_shell(j+1:j+2, 1);
                dataFoldD(count).y = blue_shell(j+1:j+2, 2);
                dataFoldD(count).color = blue;

            end
        else
            for j = (n-4):2:2*(n-4)-2
                
                % Indexing for horizontal offset
                index = (n/2) + ((j-(n-4))/2) + 1;
                
                % Lengths index
                lengthsindex = (n/2) + ((j-(n-4))/2) + 2;
                
                % Populate array
                blue_shell(j+1, 1) = index*ls - ls/2;
                blue_shell(j+1, 2) = h1 + max(lengths);
                blue_shell(j+2, 1) = index*ls;
                blue_shell(j+2, 2) = h1 + 2*max(lengths) - lengths(lengthsindex);
                
                % Increase count, add to structure, and plot
                count = count + 1;
                dataFoldD(count).x = blue_shell(j+1:j+2, 1);
                dataFoldD(count).y = blue_shell(j+1:j+2, 2);
                dataFoldD(count).color = blue;

            end
        end
    end
end


% Create "weave" pattern (pt.7) - top shell (orange)
% ------------------------------------------------------------------

% Nearly identical to process from before
if n > 4
    % Initialization
    orange_shell = zeros(2*(n-4), 2);
    
    for i = 1:2 % Identifies first or second section
        if i == 1
            for j = 0:2:(n-4)-2
                
                % Indexing for horizontal offset
                index = (j/2) + 1;
                
                % Lengths index
                lengthsindex = (j/2) + 2;
                
                % Populate array
                orange_shell(j+1, 1) = index*ls - ls/2;
                orange_shell(j+1, 2) = h1 + max(lengths);
                orange_shell(j+2, 1) = index*ls;
                orange_shell(j+2, 2) = h1 + 2*max(lengths) - lengths(lengthsindex);
                
                % Increase count, add to structure, and plot
                count = count + 1;
                dataFoldD(count).x = orange_shell(j+1:j+2, 1);
                dataFoldD(count).y = orange_shell(j+1:j+2, 2);
                dataFoldD(count).color = orange;
                
            end
        else
            for j = (n-4):2:2*(n-4)-2
                
                % Indexing for horizontal offset
                index = (n/2) + ((j-(n-4))/2) + 1;
                
                % Lengths index
                lengthsindex = (n/2) + ((j-(n-4))/2) + 2;
                
                % Populate array
                orange_shell(j+1, 1) = index*ls;
                orange_shell(j+1, 2) = h1 + 2*max(lengths) - lengths(lengthsindex);
                orange_shell(j+2, 1) = index*ls + ls/2;
                orange_shell(j+2, 2) = h1 + max(lengths);
                
                % Increase count, add to structure, and plot
                count = count + 1;
                dataFoldD(count).x = orange_shell(j+1:j+2, 1);
                dataFoldD(count).y = orange_shell(j+1:j+2, 2);
                dataFoldD(count).color = orange;

            end
        end
    end
end

% Recursive Sink Gadget Addition (For n = 4)
% -----------------------------------------------------------------

if nz > 1 && n == 4
    
    % Begin by determining the slope of the diagonal lines. 
    lineslope = 2*max(lengths) / ls;
    
    % We know that the size of the recursion array will simply be 10x the
    % value of nz (since it is simply a rectangle in each occurrence)
    recursion = zeros(10*(nz-1), 2);
    
    % Determine points on x pattern where box outline will fall based on
    % the value of nz 
    
    % First half
    for i = 1:(nz-1)
        
        % Increase count initially
        count = count + 1;
        
        % Index for iteration
        index = 5*(i-1)+1;
        
        % Populate the values of recursion based on the value of lineslope
        % and the respective position on the x axis
        recursion(index, 1) = (i/nz)*(ls/2);
        recursion(index, 2) = h1 + 2*max(lengths) - (lineslope*recursion(index, 1));
        
        recursion(index+1, 1) = (i/nz)*(ls/2);
        recursion(index+1, 2) = h1 + lineslope*recursion(index+1, 1);
        
        recursion(index+2, 1) = ls - (i/nz)*(ls/2);
        recursion(index+2, 2) = h1 + 2*max(lengths) - lineslope*recursion(index+2, 1);
        
        recursion(index+3, 1) = ls - (i/nz)*(ls/2);
        recursion(index+3, 2) = h1 + lineslope*recursion(index+3, 1);
        
        recursion(index+4, 1) = recursion(index, 1);
        recursion(index+4, 2) = recursion(index, 2);
        
        % Storing information to dataFoldD
        dataFoldD(count).x = recursion(index:index+4, 1);
        dataFoldD(count).y = recursion(index:index+4, 2);
        dataFoldD(count).color = blue;
        
    end
    
    % Second half
    for i = 1:(nz-1)

        % Increase count initially
        count = count + 1;

        % Index for iteration
        index = 5*(i-1)+(5*(nz-1))+1;
        
        % Determine vertical offset from slope and horizontal offset
        vert_offset = lineslope*2*ls;

        % Populate the values of recursion based on the value of lineslope
        % and the respective position on the x axis
        recursion(index, 1) = (2*ls) + (i/nz)*(ls/2);
        recursion(index, 2) = vert_offset + h1 + 2*max(lengths) - (lineslope*recursion(index, 1));

        recursion(index+1, 1) = (2*ls) + (i/nz)*(ls/2);
        recursion(index+1, 2) = h1 - vert_offset + lineslope*recursion(index+1, 1);

        recursion(index+2, 1) = (3*ls) - (i/nz)*(ls/2);
        recursion(index+2, 2) = vert_offset + h1 + 2*max(lengths) - lineslope*recursion(index+2, 1);

        recursion(index+3, 1) = (3*ls) - (i/nz)*(ls/2);
        recursion(index+3, 2) = h1 - vert_offset + lineslope*recursion(index+3, 1);

        recursion(index+4, 1) = recursion(index, 1);
        recursion(index+4, 2) = recursion(index, 2);

        % Storing information to dataFoldD
        dataFoldD(count).x = recursion(index:index+4, 1);
        dataFoldD(count).y = recursion(index:index+4, 2);
        dataFoldD(count).color = blue;

    end    
    
end

% Recursive Sink Gadget Addition (For all greater than n = 4)
% -----------------------------------------------------------------

if nz > 1 && n >= 6

% --------------------- General Initialization -------------------------
    
    % First, determine the slope of the diagonal lines. This value is
    % important to know, as many of the points which will be added are
    % located along these lines
    lineslope = 2*max(lengths) / ls;
    
    % We also need to determine the diagonal slopes for the diamond
    % patterns. We can store these in the diamondslope array.
    diamondslope = zeros((n-4)*2, 1);
    
    % Determine vertical offset from lineslope and horizontal offset
    vert_offset = lineslope*(n/2)*ls;
    
    % Populate diamondslope based on the values contained in slopes
    
    % First half of pattern
    for i = 0:2:(n-4)-2
        
        index = (i/2) + 2;
        
        % Slopes for the first half of array
        diamondslope(i+1) = (lengths(index) - max(lengths)) / (ls/2); 
        diamondslope(i+2) = (max(lengths) - lengths(index)) / (ls/2);
        
    end
    
    % Second half of pattern
     for i = (n-4):2:((n-4)*2)-2
         
         index = (i-(n-4))/2 + 2;
         
         % Slopes for the second half of array
         diamondslope(i+1) = (lengths(index) - max(lengths)) / (ls/2);
         diamondslope(i+2) = (max(lengths) - lengths(index)) / (ls/2);        
        
    end
    
    % Initialize array to hold plot values
    recursion = zeros(((nz-1)*9)*(n-2), 2);

% ----------------------- Left "Endcaps" -------------------------------
    
    % Begin by defining and plotting "endcap" zones. Similar process to
    % situation for n = 4. For cases of n = 6, this will essentially
    % populate the entire recursive joint
    for i = 1:(nz-1)
        
        % Determine shift value, which will be used in placing
        % diamond-intersecting points.

        % Use any point on line to determine angle based on diamondslope
        v1 = [ls, max(lengths)+h1, 0] - [ls/2, max(lengths)+h1, 0];
        v2 = [ls, (diamondslope(1)*ls)+h1+max(lengths), 0] - [ls/2, max(lengths)+h1, 0];
        shiftangle = atan2(norm(cross(v1, v2)), dot(v1, v2));
                        
        % Determine shift value (x value in supp. diagram)
        shift = (((nz-i)/nz)*(ls/2))*(1/cos(shiftangle));
        
        % Use lineslope and diamondslope values to determine h, which is
        % the perpendicular from the midline to the inclined inner revolute
        % point
        h = ((lineslope)*shift) / (1 + (lineslope*diamondslope(2)));
        
        % Increase count initially
        count = count + 1;
        
        % Left segments
        index1 = 9*(i-1) + 1;
        index2 = 9*(i-1) + 1 + ((nz-1)*9)*((n/2)-1);
        
        % Populate values of recursion based on index1 and index2
        
        % index1
        recursion(index1, 1) = ls/2 + (h/lineslope);
        recursion(index1, 2) = h1 + (lineslope*recursion(index1, 1));
        
        % Double points included for facilitation of consistency of array
        % size (helps with looping)
        recursion(index1+1, 1) = (i/nz)*(ls/2);
        recursion(index1+1, 2) = h1 + 2*max(lengths) - (lineslope*recursion(index1+1, 1));
        
        recursion(index1+2, 1) = (i/nz)*(ls/2);
        recursion(index1+2, 2) = h1 + 2*max(lengths) - (lineslope*recursion(index1+2, 1));
        
        recursion(index1+3, 1) = (i/nz)*(ls/2);
        recursion(index1+3, 2) = h1 + (lineslope*recursion(index1+3, 1));
        
        recursion(index1+4, 1) = (i/nz)*(ls/2);
        recursion(index1+4, 2) = h1 + (lineslope*recursion(index1+4, 1));
        
        recursion(index1+5, 1) = ls/2 + (h/lineslope);
        recursion(index1+5, 2) = h1 + 2*max(lengths) - (lineslope*recursion(index1+5, 1));
        
        % For points 7 and 8, determine if the diamondslope value will
        % correspond to a point that is less than the value of ls.
        % otherwise, the value sets to ls.
        
        % Determine slope perp. to diamond region
        ascslope = -1/diamondslope(1);
        
        % Determine horizontal intersection at ls
        ascoffset = (recursion(index1+5, 2)) - recursion(index1+5, 1)*ascslope;
        ascintersect = ((max(lengths) + h1)  - ascoffset) / ascslope;
        
        if ascintersect < ls
            
            % If the intersection occurs before frame limit, graph normally
            % without factoring in cutoff.
            recursion(index1+6, 1) = ascintersect;
            recursion(index1+6, 2) = max(lengths) + h1;
            
            recursion(index1+7, 1) = ascintersect;
            recursion(index1+7, 2) = max(lengths) + h1;            
            
        elseif ascintersect >= ls
            
            % If the intersection occurs after frame limit, graph only up
            % to the frame limit, then graph vertical line. 
            recursion(index1+6, 1) = ls;
            recursion(index1+6, 2) = ascslope*ls + ascoffset;
            
            % Find offset from bottom line
            differential = recursion(index1+6, 2) - recursion(index1+5, 2);
            
            recursion(index1+7, 1) = ls;
            recursion(index1+7, 2) = recursion(index1, 2) - differential;
            
        end
        
        recursion(index1+8, 1) = recursion(index1, 1);
        recursion(index1+8, 2) = recursion(index1, 2);
        
        % Storing information to dataFoldD
        dataFoldD(count).x = recursion(index1:index1+8, 1);
        dataFoldD(count).y = recursion(index1:index1+8, 2);
        dataFoldD(count).color = blue;
        
        % Increase count
        count = count + 1;
        
        % For left endcap of right pattern (index2)
        recursion(index2, 1) = ls*((n/2)+1) - (ls/2) + (h/lineslope);
        recursion(index2, 2) = -vert_offset + h1 + (lineslope*recursion(index2, 1));
        
        recursion(index2+1, 1) = (i/nz)*(ls/2) + ls*(n/2);
        recursion(index2+1, 2) = vert_offset + h1 + 2*max(lengths) - (lineslope*recursion(index2+1, 1));
        
        recursion(index2+2, 1) = (i/nz)*(ls/2) + ls*(n/2);
        recursion(index2+2, 2) = vert_offset + h1 + 2*max(lengths) - (lineslope*recursion(index2+2, 1));
        
        recursion(index2+3, 1) = (i/nz)*(ls/2) + ls*(n/2);
        recursion(index2+3, 2) = -vert_offset + h1 + (lineslope*recursion(index2+3, 1));
        
        recursion(index2+4, 1) = (i/nz)*(ls/2) + ls*(n/2);
        recursion(index2+4, 2) = -vert_offset + h1 + (lineslope*recursion(index2+4, 1));
        
        recursion(index2+5, 1) = ls*((n/2)+1) - (ls/2) + (h/lineslope);
        recursion(index2+5, 2) = vert_offset + h1 + 2*max(lengths) - (lineslope*recursion(index2+5, 1));
        
        % For points 7 and 8, determine if the diamondslope value will
        % correspond to a point that is less than the value of ls.
        % otherwise, the value sets to ls.
        
        % Determine slope perp. to diamond region
        ascslope = -1/diamondslope(n-3);
        
        % Determine horizontal intersection at ls
        ascoffset = (recursion(index2+5, 2)) - recursion(index2+5, 1)*ascslope;
        ascintersect = ((max(lengths) + h1)  - ascoffset) / ascslope;
        
        if ascintersect < ((n/2)+1)*ls
            
            % If the intersection occurs before frame limit, graph normally
            % without factoring in cutoff.
            recursion(index2+6, 1) = ascintersect;
            recursion(index2+6, 2) = max(lengths) + h1;
            
            recursion(index2+7, 1) = ascintersect;
            recursion(index2+7, 2) = max(lengths) + h1;            
            
        elseif ascintersect >= ((n/2)+1)*ls
            
            % If the intersection occurs after frame limit, graph only up
            % to the frame limit, then graph vertical line. 
            recursion(index2+6, 1) = ((n/2)+1)*ls;
            recursion(index2+6, 2) = ascslope*recursion(index2+6, 1) + ascoffset;
            
            % Find offset from bottom line
            differential = recursion(index2+6, 2) - recursion(index2+5, 2);
            
            recursion(index2+7, 1) = ((n/2)+1)*ls;
            recursion(index2+7, 2) = recursion(index2, 2) - differential;
            
        end
        
        recursion(index2+8, 1) = recursion(index2, 1);
        recursion(index2+8, 2) = recursion(index2, 2);
        
        % Storing information to dataFoldD
        dataFoldD(count).x = recursion(index2:index2+8, 1);
        dataFoldD(count).y = recursion(index2:index2+8, 2);
        dataFoldD(count).color = blue; 
        
% ------------------------- Right "Endcaps" ------------------------------
        
        % Increase count
        count = count + 1;
       
        % Indexing for right endcaps 
        index3 = 9*(i-1) + 9*(nz-1)*((n/2)-2) + 1;
        index4 = 9*(i-1) + 9*(nz-1)*((n-3)) + 1;
        
        % Populate values of recursion based on index3 and index4
        
        % The y intersect in the positive regime (for neg. slope)
        index3offsetpos = (h1) + lineslope*((n/2)-1)*ls;
        
        % The y intersect in the negative regime (for pos. slope)
        index3offsetneg = (h1 + 2*max(lengths)) - lineslope*((n/2)-1)*ls;
        
        % index3
        recursion(index3, 1) = ls/2 - (h/lineslope) + ((n/2)-2)*ls;
        recursion(index3, 2) = index3offsetpos - (lineslope*recursion(index3, 1));
        
        recursion(index3+1, 1) = ((n/2)-1)*ls - (i/nz)*(ls/2);
        recursion(index3+1, 2) = index3offsetneg + (lineslope*recursion(index3+1, 1));
        
        recursion(index3+2, 1) = ((n/2)-1)*ls - (i/nz)*(ls/2);
        recursion(index3+2, 2) = index3offsetneg + (lineslope*recursion(index3+2, 1));
        
        recursion(index3+3, 1) = ((n/2)-1)*ls - (i/nz)*(ls/2);
        recursion(index3+3, 2) = index3offsetpos - (lineslope*recursion(index3+3, 1));
        
        recursion(index3+4, 1) = ((n/2)-1)*ls - (i/nz)*(ls/2);
        recursion(index3+4, 2) = index3offsetpos - (lineslope*recursion(index3+4, 1));
        
        recursion(index3+5, 1) = ls/2 - (h/lineslope) + ((n/2)-2)*ls;
        recursion(index3+5, 2) = index3offsetneg + (lineslope*recursion(index3+5, 1));
               
        % For points 7 and 8, determine if the value dictated by
        % diamondslope will cross over the left-side ls border. 
        
        % Determine slope perp. to diamond region
        descslope = -1/diamondslope(n-4);
        
        % Determine horizontal intersection at ls
        descoffset = (recursion(index3+5, 2)) - recursion(index3+5, 1)*descslope;
        descintersect = ((max(lengths) + h1)  - descoffset) / descslope;
        
        if descintersect > ((n/2)-2)*ls
            
            % If the intersection occurs before frame limit, graph normally
            % without factoring in cutoff.
            recursion(index3+6, 1) = descintersect;
            recursion(index3+6, 2) = max(lengths) + h1;
            
            recursion(index3+7, 1) = descintersect;
            recursion(index3+7, 2) = max(lengths) + h1;            
            
        elseif descintersect <= ((n/2)-2)*ls
            
            % If the intersection occurs after frame limit, graph only up
            % to the frame limit, then graph vertical line. 
            recursion(index3+6, 1) = ((n/2)-2)*ls;
            recursion(index3+6, 2) = descslope*recursion(index3+6, 1) + descoffset;
            
            % Find offset from bottom line
            differential = recursion(index3+6, 2) - recursion(index3+5, 2);
            
            recursion(index3+7, 1) = ((n/2)-2)*ls;
            recursion(index3+7, 2) = recursion(index3, 2) - differential;
            
        end
        
        recursion(index3+8, 1) = recursion(index3, 1);
        recursion(index3+8, 2) = recursion(index3, 2);
        
        % Storing information to dataFoldD
        dataFoldD(count).x = recursion(index3:index3+8, 1);
        dataFoldD(count).y = recursion(index3:index3+8, 2);
        dataFoldD(count).color = blue;
        
        % Increase count
        count = count + 1;
        
        % The y intersect in the positive regime (for neg. slope)
        index4offsetpos = (h1) + lineslope*(n-1)*ls;
        
        % The y intersect in the negative regime (for pos. slope)
        index4offsetneg = (h1 + 2*max(lengths)) - lineslope*(n-1)*ls;
        
        % index3
        recursion(index4, 1) = ls/2 - (h/lineslope) + (n-2)*ls;
        recursion(index4, 2) = index4offsetpos - (lineslope*recursion(index4, 1));
        
        recursion(index4+1, 1) = (n-1)*ls - (i/nz)*(ls/2);
        recursion(index4+1, 2) = index4offsetneg + (lineslope*recursion(index4+1, 1));
        
        recursion(index4+2, 1) = (n-1)*ls - (i/nz)*(ls/2);
        recursion(index4+2, 2) = index4offsetneg + (lineslope*recursion(index4+2, 1));
        
        recursion(index4+3, 1) = (n-1)*ls - (i/nz)*(ls/2);
        recursion(index4+3, 2) = index4offsetpos - (lineslope*recursion(index4+3, 1));
        
        recursion(index4+4, 1) = (n-1)*ls - (i/nz)*(ls/2);
        recursion(index4+4, 2) = index4offsetpos - (lineslope*recursion(index4+4, 1));
        
        recursion(index4+5, 1) = ls/2 - (h/lineslope) + (n-2)*ls;
        recursion(index4+5, 2) = index4offsetneg + (lineslope*recursion(index4+5, 1));
               
        % For points 7 and 8, determine if the value dictated by
        % diamondslope will cross over the left-side ls border. 
        
        % Determine slope perp. to diamond region
        descslope = -1/diamondslope(end);
        
        % Determine horizontal intersection at ls
        descoffset = (recursion(index4+5, 2)) - recursion(index4+5, 1)*descslope;
        descintersect = ((max(lengths) + h1)  - descoffset) / descslope;
        
        if descintersect > (n-2)*ls
            
            % If the intersection occurs before frame limit, graph normally
            % without factoring in cutoff.
            recursion(index4+6, 1) = descintersect;
            recursion(index4+6, 2) = max(lengths) + h1;
            
            recursion(index4+7, 1) = descintersect;
            recursion(index4+7, 2) = max(lengths) + h1;            
            
        elseif descintersect <= (n-2)*ls
            
            % If the intersection occurs after frame limit, graph only up
            % to the frame limit, then graph vertical line. 
            recursion(index4+6, 1) = (n-2)*ls;
            recursion(index4+6, 2) = descslope*recursion(index4+6, 1) + descoffset;
            
            % Find offset from bottom line
            differential = recursion(index3+6, 2) - recursion(index3+5, 2);
            
            recursion(index4+7, 1) = (n-2)*ls;
            recursion(index4+7, 2) = recursion(index4, 2) - differential;
            
        end
        
        recursion(index4+8, 1) = recursion(index4, 1);
        recursion(index4+8, 2) = recursion(index4, 2);
        
        % Storing information to dataFoldD
        dataFoldD(count).x = recursion(index4:index4+8, 1);
        dataFoldD(count).y = recursion(index4:index4+8, 2);
        dataFoldD(count).color = blue;
      
    end   
    
end

% ----------------- Provisional Additions for n > 6 ---------------------

if nz > 1 && n >= 8
    
    % Determine number of ls segments to iterate through
    num = n - 6;
    
    % Loop through to assign the layered recursion values (left sect.)
    for i = 1:(nz-1)
        
        % ---------------- Left Half of Pattern ------------------------
        
        for j = 1:num/2 % Iteration through first half of pattern
            
            % New indexing for diamonds (as they scale by 2)
            diamondindexeven = 2*j;
            diamondindexodd = (2*j) + 1;
            
            % Determine shift value, which will be used in placing
            % diamond-intersecting points.

            % Use any point on line to determine angles based on
            % diamondslope values
            v11 = [ls, max(lengths)+h1, 0] - [ls/2, max(lengths)+h1, 0];
            v21 = [ls, (diamondslope(diamondindexeven-1)*ls)+h1+max(lengths), 0] - [ls/2, max(lengths)+h1, 0];
            shiftangle1 = atan2(norm(cross(v11, v21)), dot(v11, v21));
            
            v12 = [ls, max(lengths)+h1, 0] - [ls/2, max(lengths)+h1, 0];
            v22 = [ls, (diamondslope(diamondindexodd)*ls)+h1+max(lengths),0] - [ls/2, max(lengths)+h1, 0];
            shiftangle2 = atan2(norm(cross(v12, v22)), dot(v12, v22));

            % Determine shift values (x value in supp. diagram)
            shift1 = (((nz-i)/nz)*(ls/2))*(1/cos(shiftangle1));
            shift2 = (((nz-i)/nz)*(ls/2))*(1/cos(shiftangle2));
            
            % Use lineslope and diamondslope values to determine h, which is
            % the perpendicular from the midline to the inclined inner revolute
            % point
            hshift1 = ((lineslope)*shift1) / (1 + (lineslope*diamondslope(diamondindexeven)));
            hshift2 = ((lineslope)*shift2) / (1 + (lineslope*-diamondslope(diamondindexodd)));
                        
            index = (j-1)*9 + (nz-1)*9;
            
            % Determine y offsets (b values)
            pos = (h1) + lineslope*(j+1)*ls;
        
            % The y intersect in the negative regime (for pos. slope)
            neg = (h1 + 2*max(lengths)) - lineslope*(j+1)*ls;
            
            
            % Populate points
            recursion(index, 1) = ls/2 + (hshift2/lineslope) + j*ls;
            recursion(index, 2) = neg + lineslope*recursion(index, 1);
            
            recursion(index+1, 1) = ls/2 - (hshift1/lineslope) + j*ls;
            recursion(index+1, 2) = pos - lineslope*recursion(index+1, 1);
            
            recursion(index+2, 1) = ls/2 - shift1 + j*ls;
            recursion(index+2, 2) = max(lengths) + h1;
            
            recursion(index+3, 1) = ls/2 - shift1 + j*ls;
            recursion(index+3, 2) = max(lengths) + h1;
            
            recursion(index+4, 1) = ls/2 - (hshift1/lineslope) + j*ls;
            recursion(index+4, 2) = neg + lineslope*recursion(index+4, 1);
            
            recursion(index+5, 1) = ls/2 + (hshift2/lineslope) + j*ls;
            recursion(index+5, 2) = pos - lineslope*recursion(index+5, 1);
            
            recursion(index+6, 1) = ls/2 + shift2 + j*ls;
            recursion(index+6, 2) = max(lengths) + h1;
            
            recursion(index+7, 1) = ls/2 + shift2 + j*ls;
            recursion(index+7, 2) = max(lengths) + h1;
            
            recursion(index+8, 1) = recursion(index, 1);
            recursion(index+8, 2) = recursion(index, 2);
                            
            % Storing and Plotting
            count = count + 1;
            dataFoldD(count).x = recursion(index:index+8, 1);
            dataFoldD(count).y = recursion(index:index+8, 2);
            dataFoldD(count).color = blue;
            
        end 
        
        % ------------------ Right Half of Pattern ---------------------
        
        for j = 1:num/2 % Iteration through second half of pattern
            
            % New indexing for diamonds (as they scale by 2)
            diamondindexeven = 2*j + (n-4);
            diamondindexodd = 2*j + (n-4) + 1;
            
            % Determine shift values, which will be used in placing
            % diamond-intersecting points.

            % Use any point on line to determine angles based on
            % diamondslope values
            v11 = [ls, max(lengths)+h1, 0] - [ls/2, max(lengths)+h1, 0];
            v21 = [ls, (diamondslope(diamondindexeven-1)*ls)+h1+max(lengths), 0] - [ls/2, max(lengths)+h1, 0];
            shiftangle1 = atan2(norm(cross(v11, v21)), dot(v11, v21));
            
            v12 = [ls, max(lengths)+h1, 0] - [ls/2, max(lengths)+h1, 0];
            v22 = [ls, (diamondslope(diamondindexodd)*ls)+h1+max(lengths),0] - [ls/2, max(lengths)+h1, 0];
            shiftangle2 = atan2(norm(cross(v12, v22)), dot(v12, v22));

            % Determine shift values (x value in supp. diagram)
            shift1 = (((nz-i)/nz)*(ls/2))*(1/cos(shiftangle1));
            shift2 = (((nz-i)/nz)*(ls/2))*(1/cos(shiftangle2));
            
            % Use lineslope and diamondslope values to determine h, which is
            % the perpendicular from the midline to the inclined inner revolute
            % point
            hshift1 = ((lineslope)*shift1) / (1 + (lineslope*diamondslope(diamondindexeven)));
            hshift2 = ((lineslope)*shift2) / (1 + (lineslope*-diamondslope(diamondindexodd)));
                        
            index = (j-1)*9 + (n/2)*9*(nz-1) + 1;
            
            % Determine ls offset value
            offset = (n/2) + 1;
            
            % Determine y offsets (b values)
            pos = (h1) + lineslope*(j+offset)*ls;
        
            % The y intersect in the negative regime (for pos. slope)
            neg = (h1 + 2*max(lengths)) - lineslope*(j+offset)*ls;
            
            % Populate "corner" points first (on lineslope x)
            recursion(index, 1) = ls/2 + (hshift2/lineslope) + (j+offset-1)*ls;
            recursion(index, 2) = neg + lineslope*recursion(index, 1);
            
            recursion(index+1, 1) = ls/2 - (hshift1/lineslope) + (j+offset-1)*ls;
            recursion(index+1, 2) = pos - lineslope*recursion(index+1, 1);
            
            recursion(index+2, 1) = ls/2 - shift1 + (j+offset-1)*ls;
            recursion(index+2, 2) = max(lengths) + h1;
            
            recursion(index+3, 1) = ls/2 - shift1 + (j+offset-1)*ls;
            recursion(index+3, 2) = max(lengths) + h1;
            
            recursion(index+4, 1) = ls/2 - (hshift1/lineslope) + (j+offset-1)*ls;
            recursion(index+4, 2) = neg + lineslope*recursion(index+4, 1);
            
            recursion(index+5, 1) = ls/2 + (hshift2/lineslope) + (j+offset-1)*ls;
            recursion(index+5, 2) = pos - lineslope*recursion(index+5, 1);
            
            recursion(index+6, 1) = ls/2 + shift2 + (j+offset-1)*ls;
            recursion(index+6, 2) = max(lengths) + h1;
            
            recursion(index+7, 1) = ls/2 + shift2 + (j+offset-1)*ls;
            recursion(index+7, 2) = max(lengths) + h1;
            
            recursion(index+8, 1) = recursion(index, 1);
            recursion(index+8, 2) = recursion(index, 2);
                       
            % Storing and Plotting
            count = count + 1;
            dataFoldD(count).x = recursion(index:index+8, 1);
            dataFoldD(count).y = recursion(index:index+8, 2);
            dataFoldD(count).color = blue;
            
        end 
                
    end
  
end

% Add section for overlap fold
% -------------------------------------------------------------------

% % Orange Shell
% 
% if n > 4
%     
%     % Increase count
%     count = count + 1;
% 
%     % Initialize 2x2 array
%     overlap_shell = zeros(2, 2);
%     overlap_shell(1, 1) = orange_shell(1, 1) + n*ls;
%     overlap_shell(1, 2) = orange_shell(1, 2);
%     overlap_shell(2, 1) = orange_shell(2, 1) + n*ls;
%     overlap_shell(2, 2) = orange_shell(2, 2);
% 
%     % Store to data structure 
%     dataFoldD(count).x = overlap_shell(1:2, 1);
%     dataFoldD(count).y = overlap_shell(1:2, 2);
%     dataFoldD(count).color = orange;
% 
%     % Plot
%     plot(dataFoldD(count).x, dataFoldD(count).y, 'color', ...
%         dataFoldD(count).color)
% 
% end
% 
% % Vertical Blue Line
% 
% % Increase count
% count = count + 1;
% 
% % Initialize 2x2 array
% overlap_blue = zeros(2, 2);
% overlap_blue(1, 1) = n*ls;
% overlap_blue(1, 2) = h1;
% overlap_blue(2, 1) = n*ls;
% overlap_blue(2, 2) = h1 + 2*max(lengths);
% 
% % Store to data structure 
% dataFoldD(count).x = overlap_blue(1:2, 1);
% dataFoldD(count).y = overlap_blue(1:2, 2);
% dataFoldD(count).color = blue;
% 
% % Plot
% plot(dataFoldD(count).x, dataFoldD(count).y, 'color', ...
%     dataFoldD(count).color)
% 
% % Horizontal Blue Line
% 
% % Increase count
% count = count + 1;
% 
% % Initialize 2x2 array
% overlap_horizon_blue = zeros(2, 2);
% overlap_horizon_blue(1, 1) = n*ls;
% overlap_horizon_blue(1, 2) = h1 + max(lengths);
% overlap_horizon_blue(2, 1) = n*ls + ls/2;
% overlap_horizon_blue(2, 2) = h1 + max(lengths);
% 
% % Store to data structure 
% dataFoldD(count).x = overlap_horizon_blue(1:2, 1);
% dataFoldD(count).y = overlap_horizon_blue(1:2, 2);
% dataFoldD(count).color = blue;
% 
% % Plot
% plot(dataFoldD(count).x, dataFoldD(count).y, 'color', ...
%     dataFoldD(count).color)
% 
% % Horizontal Orange Line
% 
% % Increase count
% count = count + 1;
% 
% % Initialize 2x2 array
% overlap_horizon_orange = zeros(2, 2);
% overlap_horizon_orange(1, 1) = n*ls + ls/2;
% overlap_horizon_orange(1, 2) = h1 + max(lengths);
% overlap_horizon_orange(2, 1) = (n+1)*ls;
% overlap_horizon_orange(2, 2) = h1 + max(lengths);
% 
% % Store to data structure 
% dataFoldD(count).x = overlap_horizon_orange(1:2, 1);
% dataFoldD(count).y = overlap_horizon_orange(1:2, 2);
% 
% if n > 4
%     dataFoldD(count).color = orange;
% else
%     dataFoldD(count).color = blue;
% end
% 
% 
% % Plot
% plot(dataFoldD(count).x, dataFoldD(count).y, 'color', ...
%     dataFoldD(count).color)

% ------------------ Yellow Lines for 3D Rendering ----------------------

% % Left Panel
% xcoordsL = [((n/2)-1)*ls, (n/2)*ls, ((n/2)-1)*ls];
% ycoordsL = [0, h1 + lmax, h1 + h2 + 2*lmax];
% 
% % Increase count, store to data structure, and plot
% count = count + 1;
% dataFoldD(count).x = xcoordsL;
% dataFoldD(count).y = ycoordsL;
% dataFoldD(count).color = yellow;
% 
% plot(dataFoldD(count).x, dataFoldD(count).y, 'color', ...
%     dataFoldD(count).color)
% 
% 
% % Right Panel
% xcoordsR = [(n-1)*ls, n*ls, (n-1)*ls];
% ycoordsR = [0, h1 + lmax, h1 + h2 + 2*lmax];
% 
% % Increase count, store to data structure, and plot
% count = count + 1;
% dataFoldD(count).x = xcoordsR;
% dataFoldD(count).y = ycoordsR;
% dataFoldD(count).color = yellow;
% 
% plot(dataFoldD(count).x, dataFoldD(count).y, 'color', ...
%     dataFoldD(count).color)

% Finally, plot labeling and adjustments
% ------------------------------------------------------------------

% Label the plot for clarity
% title({
%     ('Origami Schematic 2.A for Provided Parameters:')
%     ['[r = ' num2str(r) ', n = ' num2str(n) ', theta = ', num2str(theta_m) ']']
%     })

% daspect([1 1 1])

m = 0;
lmax = h1 + 2*lmax + h2;

end