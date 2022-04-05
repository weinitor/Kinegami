function [msum, lmax_sum, infostruct, Struct1, Struct2] =  ...
    DataFoldDuplicate(a, Struct1, Struct2, infostruct, index, msum, lmax_sum, triple)
% DATAFOLDUPLICATE - please describe this function.
%   infostruct - a data structure that includes all the information needed
%                for the construction of the full schematic.
%   index      - allows us to know what section of the template we are on. 
%                Must be an odd number.
%   msum       - a cumulative counter used for resetting along the x axis 
%                that carries over for each set of appended structures.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last Edited 6/22/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.txt for detail.


% First, check to ensure that fold type is not "twist"
if strcmp(infostruct(index).name, "Twist") == 1
    name = 1;
else
    name = 0;
end

% num will be used to prevent double lines for twists that become tubes
if mod(abs(infostruct(index).m), infostruct(index).ls) < 10^-4 || ...
        mod(abs(infostruct(index).m), infostruct(index).ls) > infostruct(index).ls - 10^-4
    num = 5;
else
    num = 4;
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
            
            for i = 1:5
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end
            
            % Ignore left tube line
            for i = 6:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
            end
            
        % If appended to left (no twist)
        elseif infostruct(index).n * infostruct(index).ls > 2*pi*infostruct(index).r ...
                && name == 0
            
            % Override right tube data with left tube data, then adjust
            % according to specified offset
            Struct1(6).x = Struct1(5).x;
            Struct1(6).y = Struct1(5).y;
            
            for i = 1:5
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end
            
            for i = 6:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
            end
        
        % If appended to right (twist)
        elseif infostruct(index).n * infostruct(index).ls <= 2*pi*infostruct(index).r ...
                && name == 1
            
            for i = 1:num
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end
            
            for i = num+1:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
            end
        
        % If appended to left (twist)
        else
            
            for i = 1:num
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end
            
            for i = num+1:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
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
            
            for i = 1:5
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end

            % Ignore left tube line
            for i = 6:size(Struct1, 2)

                % Factor in offset for every location along the structure
                Struct1(i).x = Struct1(i).x + msum + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                
                % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)

            end

        % If appended to left (not twist)
        elseif infostruct(index).n * infostruct(index).ls + msum > 2*pi*infostruct(index).r ...
                && name == 0

            % Reassign left tube data to right tube index and specify offsets
             Struct1(6).x = Struct1(5).x;
             Struct1(6).y = Struct1(5).y; 
             
            for i = 1:5
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end

             % Ignore the right tube line
             for i = 6:size(Struct1, 2)

                 % Factor in offset for every location along the structure
                 Struct1(i).x = Struct1(i).x + msum - (infostruct(index).n * infostruct(index).ls);
                 Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                 
                 % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)

             end
        
        % If appended to right (twist)
        elseif infostruct(index).n * infostruct(index).ls + msum <= 2*pi*infostruct(index).r ...
                && name == 1
            
            for i = 1:num
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end

            for i = num+1:size(Struct1, 2)

                % Factor in offset for every location along the structure
                Struct1(i).x = Struct1(i).x + msum + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                
                % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)

            end
         
        % If appended to left (twist)
        else
            
            for i = 1:num
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end
            
            for i = num+1:size(Struct1, 2)

                % Factor in offset for every location along the structure
                Struct1(i).x = Struct1(i).x + msum - (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                
                % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)

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
            
            for i = 1:5
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end
            
            % Ignore left tube line
            for i = 6:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)       
                
            end
            
            % Reassign left tube data to right tube index and specify offsets
            Struct2(6).x = Struct2(5).x;
            Struct2(6).y = Struct2(5).y;
            
            for i = 1:5
                
                % Assign null values for purpose of DXF Generation
                Struct2(i).x = [];
                Struct2(i).y = [];
                
            end
            
            % Plot left component of triplet
            for i = 6:size(Struct2, 2)
                
                Struct2(i).x = Struct2(i).x - (infostruct(index).n * infostruct(index).ls);
                Struct2(i).y = Struct2(i).y;
                
                % Plotting
                plot(a, Struct2(i).x, Struct2(i).y, 'color', Struct2(i).color)
                
            end
            
        end
        
        % If appended to right (twist)
        if name == 1
            
            for i = 1:num
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end
            
            % Ignore boundary lines
            for i = num+1:size(Struct1, 2)
                
                Struct1(i).x = Struct1(i).x + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y;
                
                % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
                
            end
            
            for i = 1:num
                
                % Assign null values for purpose of DXF Generation
                Struct2(i).x = [];
                Struct2(i).y = [];
                
            end
            
            % Left triplet component
            for i = num+1:size(Struct2, 2)
                
                Struct2(i).x = Struct2(i).x - (infostruct(index).n * infostruct(index).ls);
                Struct2(i).y = Struct2(i).y;
                
                % Plot each segment with offset
                plot(a, Struct2(i).x, Struct2(i).y, 'color', Struct2(i).color)
                
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
            
            for i = 1:5
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end
            
            % Ignore left tube line
            for i = 6:size(Struct1, 2)

                % Factor in offset for every location along the structure
                Struct1(i).x = Struct1(i).x + msum + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;
                
                % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)

            end
            
            % Reassign left tube data to right tube index and specify offsets
            Struct2(6).x = Struct2(5).x;
            Struct2(6).y = Struct2(5).y;
            
           for i = 1:5
                
                % Assign null values for purpose of DXF Generation
                Struct2(i).x = [];
                Struct2(i).y = [];
                
            end

            % Left triplet component
            for i = 6:size(Struct2, 2)
                
                % Factor in offset for every location along the structure
                Struct2(i).x = Struct2(i).x + msum - (infostruct(index).n * infostruct(index).ls);
                Struct2(i).y = Struct2(i).y + infostruct(index-1).lmax + lmax_sum;
                 
                % Plot each segment as before, but with offset factored in
                plot(a, Struct2(i).x, Struct2(i).y, 'color', Struct2(i).color)
                
            end
            
        end

        % Twist
        if name == 1
            
            for i = 1:num
                
                % Assign null values for purpose of DXF Generation
                Struct1(i).x = [];
                Struct1(i).y = [];
                
            end

            % Plot on right side
            for i = num+1:size(Struct1, 2)

                % Factor in offset for every location along the structure
                Struct1(i).x = Struct1(i).x + msum + (infostruct(index).n * infostruct(index).ls);
                Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;

                % Plot each segment as before, but with offset factored in
                plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
                
            end
            
            for i = 1:num
                
                % Assign null values for purpose of DXF Generation
                Struct2(i).x = [];
                Struct2(i).y = [];
                
            end
            
            % Plot left side of triplet
            for i = num+1:size(Struct2, 2)
                
                % Factor in left offset
                Struct2(i).x = Struct2(i).x + msum - (infostruct(index).n * infostruct(index).ls);
                Struct2(i).y = Struct2(i).y + infostruct(index-1).lmax + lmax_sum;
                
                % Plot each segment with offset
                plot(a, Struct2(i).x, Struct2(i).y, 'color', Struct2(i).color)
                
            end
        end     
    
        lmax_sum = infostruct(index-1).lmax + lmax_sum;
    
    end    

end

end