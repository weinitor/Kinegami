% Duplicating Crease Pattern
% Last Edited 6/22/2021 by Lucien Peach

function [msum, lmax_sum] = DataFoldDuplicate(Struct1, infostruct, index, msum, lmax_sum, triple)
% infostruct is a data structure that includes all the information needed
% for the construction of the full schematic

% index allows us to know what section of the template we are on. Must be
% an odd number

% msum is a cumulative counter used for resetting along the x axis that
% carries over for each set of appended structures

% First, check to ensure that fold type is not "twist"
if strcmp(infostruct(index).name, "Twist") == 1
    name = 1;
else
    name = 0;
end

% Check will differ based on if index is final in array
if index == size(infostruct, 2)
    
    % Last index check
    if infostruct(index).n == infostruct(index-1).n ...
            && infostruct(index).ls == infostruct(index-1).ls
        
        alert = 0;
    else
        alert = msgbox('Section Width and Number of Vertices Must be Identical', ...
            'Error', 'error');
        
    end
    
else % Any other index check

    if infostruct(index).n == infostruct(index+1).n ...
            && infostruct(index).ls == infostruct(index+1).ls
        % Do nothing - continue with code
        alert = 0;
    else
        alert = msgbox('Section Width and Number of Vertices Must be Identical', ...
            'Error', 'error');   
    end
    
end

% A few notes: in order for the above condition to hold, the values of n
% and r must match for all tube components. Thus indexing of r does not
% matter, as long as it can stay consistent throughout loops. Same with n.
% However indexing for m and lmax does matter. 
    
if alert == 0 && strcmp(triple, 'triple') == 0
    
    % Include special case for the first block
    if index == 1

        % No mod term needed

        % Do not consider the m offset here, since it will always be 0.
        % Be sure to factor in offset from initial frame.
        
        % If appended to right (no twist)
        if infostruct(index).n * infostruct(index).ls <= 2*pi*infostruct(index).r ... 
                && name == 0            
            % Ignore left tube line
            for i = 6:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
            end
            
        % If appended to left (no twist)
        elseif infostruct(index).n * infostruct(index).ls > 2*pi*infostruct(index).r ...
                && name == 0
            
            % Override right tube data with left tube data, then adjust
            % according to specified offset
            Struct1(6).x = Struct1(5).x;
            Struct1(6).y = Struct1(5).y;
            
            for i = 6:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
            end
        
        % If appended to right (twist)
        elseif infostruct(index).n * infostruct(index).ls <= 2*pi*infostruct(index).r ...
                && name == 1
            
            for i = 5:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
            end
        
        % If appended to left (twist)
        else
            
            for i = 5:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
            end
            
        end
                
        % End of index = 1.
             
    else % If index is not equal to 1

        % Mod term
        if mod(msum + infostruct(index-1).m, 2*pi*infostruct(index).r) ... 
                < infostruct(index-1).m + msum && mod(msum + infostruct(index-1).m, ...
                2*pi*infostruct(index).r) ~= 0 

            % If the beginning of this structure is over the 2*pi*r limit,
            % adjust so that it will be reset but offset by the m offset
            % specified by the previous structure
            msum = infostruct(index-1).m;

        else
            msum = msum + infostruct(index-1).m;       
        end   

        % First fold

        % If appended to right (not twist)
        if infostruct(index).n * infostruct(index).ls + msum <= 2*pi*infostruct(index).r ...
                && name == 0

            % Ignore left tube line
            for i = 6:size(Struct1, 2)

                % Factor in offset for every location along the structure
                Struct1(i).x = Struct1(i).x + msum + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                
                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)

            end

        % If appended to left (not twist)
        elseif infostruct(index).n * infostruct(index).ls + msum > 2*pi*infostruct(index).r ...
                && name == 0

            % Reassign left tube data to right tube index and specify offsets
             Struct1(6).x = Struct1(5).x;
             Struct1(6).y = Struct1(5).y; 

             % Ignore the right tube line
             for i = 6:size(Struct1, 2)

                 % Factor in offset for every location along the structure
                 Struct1(i).x = Struct1(i).x + msum - (infostruct(index).n * infostruct(index).ls);
                 Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                 
                 % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)

             end
        
        % If appended to right (twist)
        elseif infostruct(index).n * infostruct(index).ls + msum <= 2*pi*infostruct(index).r ...
                && name == 1

            for i = 5:size(Struct1, 2)

                % Factor in offset for every location along the structure
                Struct1(i).x = Struct1(i).x + msum + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                
                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)

            end
         
        % If appended to left (twist)
        else
            
            for i = 5:size(Struct1, 2)

                % Factor in offset for every location along the structure
                Struct1(i).x = Struct1(i).x + msum - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                
                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)

            end
        end
    
    lmax_sum = infostruct(index-1).lmax + lmax_sum;
    
    end    

end

if alert == 0 && strcmp(triple, 'triple') == 1
    
    % Include special case for the first block
    if index == 1

        % No mod term needed

        % Do not consider the m offset here, since it will always be 0.
        % Be sure to factor in offset from initial frame.
        
        % No twist
        if name == 0  
            
            % Ignore left tube line
            for i = 6:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
                
                % Reset
                Struct1(i).x = Struct1(i).x - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;           
                
            end
            
            % Reassign left tube data to right tube index and specify offsets
            Struct1(6).x = Struct1(5).x;
            Struct1(6).y = Struct1(5).y;
            
            % Plot left component of triplet
            for i = 6:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plotting
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
                
            end
            
        end
        
        % If appended to right (twist)
        if name == 1
            
            % Ignore boundary lines
            for i = 5:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
                
                % Reset
                Struct1(i).x = Struct1(i).x - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
            end
            
            % Left triplet component
            for i = 5:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment with offset
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
                
            end
            
        end

        % End of index = 1.
             
    else % If index is not equal to 1

        % Mod term
        if mod(msum + infostruct(index-1).m, 2*pi*infostruct(index).r) ... 
                < infostruct(index-1).m + msum && mod(msum + infostruct(index-1).m, ...
                2*pi*infostruct(index).r) ~= 0 

            % If the beginning of this structure is over the 2*pi*r limit,
            % adjust so that it will be reset but offset by the m offset
            % specified by the previous structure
            msum = infostruct(index-1).m;

        else
            msum = msum + infostruct(index-1).m;       
        end   

        % First fold

        % If not twist
        if name == 0
            
            % Ignore left tube line
            for i = 6:size(Struct1, 2)

                % Factor in offset for every location along the structure
                Struct1(i).x = Struct1(i).x + msum + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                
                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
                
                % Reset
                Struct1(i).x = Struct1(i).x - msum - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y - infostruct(index-1).lmax - lmax_sum;

            end
            
            % Reassign left tube data to right tube index and specify offsets
            Struct1(6).x = Struct1(5).x;
            Struct1(6).y = Struct1(5).y;
            
            % Left triplet component
            for i = 6:size(Struct1, 2)
                
                % Factor in offset for every location along the structure
                Struct1(i).x = Struct1(i).x + msum - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                 
                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
                
            end
            
        end

        % Twist
        if name == 1

            % Plot on right side
            for i = 5:size(Struct1, 2)

                % Factor in offset for every location along the structure
                Struct1(i).x = Struct1(i).x + msum + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;

                % Plot each segment as before, but with offset factored in
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
                
                % Reset
                Struct1(i).x = Struct1(i).x - msum - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y - infostruct(index-1).lmax - lmax_sum;
                
            end
            
            % Plot left side of triplet
            for i = 5:size(Struct1, 2)
                
                % Factor in left offset
                Struct1(i).x = Struct1(i).x + msum - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                
                % Plot each segment with offset
                plot(Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
                
            end
        end     
    
        lmax_sum = infostruct(index-1).lmax + lmax_sum;
    
    end    

end

end