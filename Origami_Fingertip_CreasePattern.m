% Graph for crease pattern - Fingertip
% Last edited 6/11/2021 by Lucien Peach

function [FingertipFold, m, lmax] = Origami_Fingertip_CreasePattern(lengths, ls, n, h1, r, theta_m)

% Counter used for data structure indexing
count = 1;

% Identify colors
orange = [1, 0.41, 0];
blue = [0, 0, 1];
black = [0, 0, 0];

% Begin by specifying sheet boundary
% -------------------------------------------------------------------

% Determine max value of lengths array
lmax = max(lengths);

% Store boundary coordinates to array
boundarybottom = [0, 0; n*ls, 0];
boundaryleft = [0, 0; 0, h1 + 2*lmax];
boundarytop = [0, h1 + 2*lmax; n*ls, h1 + 2*lmax];
boundaryright = [n*ls, h1 + 2*lmax; n*ls, 0];

FingertipFold(count).x = boundarybottom(:, 1);
FingertipFold(count).y = boundarybottom(:, 2);
FingertipFold(count).color = black;

count = count + 1;
FingertipFold(count).x = boundaryleft(:, 1);
FingertipFold(count).y = boundaryleft(:, 2);
FingertipFold(count).color = black;

count = count + 1;
FingertipFold(count).x = boundarytop(:, 1);
FingertipFold(count).y = boundarytop(:, 2);
FingertipFold(count).color = black;

count = count + 1;
FingertipFold(count).x = boundaryright(:, 1);
FingertipFold(count).y = boundaryright(:, 2);
FingertipFold(count).color = black;

% Increase counter
count = count + 1;

FingertipFold(count).x = boundaryleft(:, 1);
FingertipFold(count).y = boundaryleft(:, 2);
FingertipFold(count).color = blue;

count = count + 1;
FingertipFold(count).x = boundaryright(:, 1);
FingertipFold(count).y = boundaryright(:, 2);
FingertipFold(count).color = blue;

% Increase counter
count = count + 1;

% Specify horizontal red lines (proximal / distal)
% ------------------------------------------------------------------

% Proximal line
proximal = [0, h1; n*ls, h1];
FingertipFold(count).x = proximal(:, 1);
FingertipFold(count).y = proximal(:, 2);
FingertipFold(count).color = blue;

% Increase counter
count = count + 1;

% Distal line
distal = [0, h1 + (2*lmax); n*ls, h1 + (2*lmax)];
FingertipFold(count).x = distal(:, 1);
FingertipFold(count).y = distal(:, 2);
FingertipFold(count).color = blue;

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
    FingertipFold(count).x = bottomtube(ii:ii+1, 1);
    FingertipFold(count).y = bottomtube(ii:ii+1, 2);
    FingertipFold(count).color = blue;

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
    toptube(jj+1, 2) = (2*lmax) + h1;
    
    % Log data to structure and add to plot. Plotting is sequential
    FingertipFold(count).x = toptube(jj:jj+1, 1);
    FingertipFold(count).y = toptube(jj:jj+1, 2);
    FingertipFold(count).color = blue;
    
end

% Cross hatch pattern
% --------------------------------------------------------------------

% The number of cross hatches will always be n-2. We can use this
% information to initialize the cross hatch array. Each cross hatch stops
% halfway, but still requires 4 points for simplicity
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
                    cross_hatch((j+2), 1) = (index-1)*ls + (ls/2);
                    
                    % Create cross hatch pattern for first batch (Y coord)
                    cross_hatch((j+1), 2) = h1;
                    cross_hatch((j+2), 2) = h1 + lmax;
                    
                    % Express in structure format and plot iteratively
                    FingertipFold(count).x = cross_hatch(j+1:j+2, 1);
                    FingertipFold(count).y = cross_hatch(j+1:j+2, 2);
                    FingertipFold(count).color = orange;
                    
                else
                    
                    % Create cross hatch pattern for second batch (X coord)
                    cross_hatch((j+3), 1) = (index-1)*ls + (ls/2);
                    cross_hatch((j+4), 1) = (index)*ls;
                    
                    % Create cross hatch pattern for second batch (Y coord)
                    cross_hatch((j+3), 2) = h1 + lmax;
                    cross_hatch((j+4), 2) = h1;
                    
                    % Express in structure format and plot iteratively
                    FingertipFold(count).x = cross_hatch(j+3:j+4, 1);
                    FingertipFold(count).y = cross_hatch(j+3:j+4, 2);
                    FingertipFold(count).color = orange;
                    
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
                    cross_hatch((j+2+offset), 1) = (index2-1)*ls + (ls/2);
                    
                    % Create cross hatch pattern for first batch (Y coord)
                    cross_hatch((j+1+offset), 2) = h1;
                    cross_hatch((j+2+offset), 2) = h1 + lmax;
                    
                    % Express in structure format and plot iteratively
                    FingertipFold(count).x = cross_hatch(j+1+offset:j+2+offset, 1);
                    FingertipFold(count).y = cross_hatch(j+1+offset:j+2+offset, 2);
                    FingertipFold(count).color = orange;
                    
                else
                    
                    % Create cross hatch pattern for second batch (X coord)
                    cross_hatch((j+3+offset), 1) = (index2-1)*ls + (ls/2);
                    cross_hatch((j+4+offset), 1) = (index2)*ls;
                    
                    % Create cross hatch pattern for second batch (Y coord)
                    cross_hatch((j+3+offset), 2) = h1 + lmax;
                    cross_hatch((j+4+offset), 2) = h1;
                    
                    % Express in structure format and plot iteratively
                    FingertipFold(count).x = cross_hatch(j+3+offset:j+4+offset, 1);
                    FingertipFold(count).y = cross_hatch(j+3+offset:j+4+offset, 2);
                    FingertipFold(count).color = orange;
                    
                end
            end            
        end 
    end
end

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

FingertipFold(count).x = weave2((o_val*4)+1:(o_val*4)+2, 1);
FingertipFold(count).y = weave2((o_val*4)+1:(o_val*4)+2, 2);
FingertipFold(count).color = orange;

% Increase count
count = count + 1;

% End segment orange horizontal, second set
weave2((o_val*4)+3, 1) = (n-1)*ls;
weave2((o_val*4)+4, 1) = n*ls;
weave2((o_val*4)+3, 2) = h1 + lengths(1);
weave2((o_val*4)+4, 2) = h1 + lengths(1);

FingertipFold(count).x = weave2((o_val*4)+3:(o_val*4)+4, 1);
FingertipFold(count).y = weave2((o_val*4)+3:(o_val*4)+4, 2);
FingertipFold(count).color = orange;

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
        FingertipFold(count).x = weave2(ii:ii+1, 1);
        FingertipFold(count).y = weave2(ii:ii+1, 2);
        FingertipFold(count).color = orange;
    
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
        FingertipFold(count).x = weave2(ii+1:ii+2, 1);
        FingertipFold(count).y = weave2(ii+1:ii+2, 2);
        FingertipFold(count).color = orange;
        
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

FingertipFold(count).x = weave1(end-7:end-6, 1);
FingertipFold(count).y = weave1(end-7:end-6, 2);
FingertipFold(count).color = blue;

% Increase count
count = count + 1;

% First line horizontal for second segment (since breaks 2x pattern)
weave1(end-3, 1) = (n/2)*ls;
weave1(end-2, 1) = (n/2)*ls + ls/2;
weave1(end-3, 2) = h1 + lengths(1);
weave1(end-2, 2) = h1 + lengths(1);

FingertipFold(count).x = weave1(end-3:end-2, 1);
FingertipFold(count).y = weave1(end-3:end-2, 2);
FingertipFold(count).color = blue;

% Increase count
count = count + 1;

% Final horizontal for segment 1
weave1(end-5, 1) = ((n/2)-1)*ls-(ls/2);
weave1(end-4, 1) = ((n/2)-1)*ls;
weave1(end-5, 2) = h1 + lengths(1);
weave1(end-4, 2) = h1 + lengths(1);

FingertipFold(count).x = weave1(end-5:end-4, 1);
FingertipFold(count).y = weave1(end-5:end-4, 2);
FingertipFold(count).color = blue;

% Increase count
count = count + 1;

% Final horizontal for second segment
weave1(end-1, 1) = (n-1)*ls-(ls/2);
weave1(end, 1) = (n-1)*ls;
weave1(end-1, 2) = h1 + lengths(1);
weave1(end, 2) = h1 + lengths(1);

FingertipFold(count).x = weave1((2*quant)-1:2*quant, 1);
FingertipFold(count).y = weave1((2*quant)-1:2*quant, 2);
FingertipFold(count).color = blue;

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

            FingertipFold(count).x = weave1(1:(quant-4), 1);
            FingertipFold(count).y = weave1(1:(quant-4), 2);
            FingertipFold(count).color = blue;

        else
            offset2 = (quant-3); % For horizontal shift

            FingertipFold(count).x = weave1(offset2:((2*quant)-8), 1);
            FingertipFold(count).y = weave1(offset2:((2*quant)-8), 2);
            FingertipFold(count).color = blue;

        end
    end
end

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

FingertipFold(count).x = vert_blue1(1:2, 1);
FingertipFold(count).y = vert_blue1(1:2, 2);
FingertipFold(count).color = blue;

% Increase count and proceed to second set (2.1 segment)
count = count + 1;

% Populate array for segment
vert_blue1(3, 1) = (n/2)*ls;
vert_blue1(3, 2) = h1;
vert_blue1(4, 1) = (n/2)*ls;
vert_blue1(4, 2) = h1 + 2*max(lengths);

FingertipFold(count).x = vert_blue1(3:4, 1);
FingertipFold(count).y = vert_blue1(3:4, 2);
FingertipFold(count).color = blue;

% Increase count and proceed to final full segment (2.2 segment)
count = count + 1;

% Populate array for segment
vert_blue1(5, 1) = (n-1)*ls;
vert_blue1(5, 2) = h1;
vert_blue1(6, 1) = (n-1)*ls;
vert_blue1(6, 2) = h1 + 2*max(lengths);

FingertipFold(count).x = vert_blue1(5:6, 1);
FingertipFold(count).y = vert_blue1(5:6, 2);
FingertipFold(count).color = blue;

% Retrace to initial segment for use within overlap patterns
count = count + 1;

% Populate array for segment
vert_blue1(7, 1) = 0;
vert_blue1(7, 2) = h1;
vert_blue1(8, 1) = 0;
vert_blue1(8, 2) = h1 + 2*max(lengths);

FingertipFold(count).x = vert_blue1(7:8, 1);
FingertipFold(count).y = vert_blue1(7:8, 2);
FingertipFold(count).color = blue;

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
                FingertipFold(count).x = vert_blue2(j+1:j+2, 1);
                FingertipFold(count).y = vert_blue2(j+1:j+2, 2);
                FingertipFold(count).color = blue;
                
                % Lengths index for varying values
                lengthsindex = (j/4) + 2;
                
                % Populate array with bottom values
                vert_blue2(j+3, 1) = index*ls;
                vert_blue2(j+3, 2) = h1;
                vert_blue2(j+4, 1) = index*ls;
                vert_blue2(j+4, 2) = lengths(lengthsindex) + h1;
                
                % Increase count and plot top values
                count = count + 1;
                FingertipFold(count).x = vert_blue2(j+3:j+4, 1);
                FingertipFold(count).y = vert_blue2(j+3:j+4, 2);
                FingertipFold(count).color = blue;

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
                FingertipFold(count).x = vert_blue2(j+1:j+2, 1);
                FingertipFold(count).y = vert_blue2(j+1:j+2, 2);
                FingertipFold(count).color = blue;
                
                % Lengths index for varying values
                lengthsindex = (n/2) + ((j-2*(n-4))/4) + 2;
                
                % Populate array with bottom values
                vert_blue2(j+3, 1) = index*ls;
                vert_blue2(j+3, 2) = h1;
                vert_blue2(j+4, 1) = index*ls;
                vert_blue2(j+4, 2) = lengths(lengthsindex) + h1;
                
                % Increase count and plot top values
                count = count + 1;
                FingertipFold(count).x = vert_blue2(j+3:j+4, 1);
                FingertipFold(count).y = vert_blue2(j+3:j+4, 2);
                FingertipFold(count).color = blue;

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
                FingertipFold(count).x = vert_orange(j+1:j+2, 1);
                FingertipFold(count).y = vert_orange(j+1:j+2, 2);
                FingertipFold(count).color = orange;

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
                FingertipFold(count).x = vert_orange(j+1:j+2, 1);
                FingertipFold(count).y = vert_orange(j+1:j+2, 2);
                FingertipFold(count).color = orange;

            end
        end
    end
end

% Finally, plot labeling and adjustments
% ------------------------------------------------------------------

% Label the plot for clarity
% title({
%     ('Origami Schematic 2.A for Provided Parameters:')
%     ['[r = ' num2str(r) ', n = ' num2str(n) ', theta = ', num2str(theta_m) ']']
%     })
% 
% daspect([1 1 1])

m = 0;
lmax = h1 + 2*lmax;

end
