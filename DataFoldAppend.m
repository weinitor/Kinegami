function [msum, lmax_sum, infostruct, Struct1] = DataFoldAppend(a, Struct1, infostruct, index, msum, lmax_sum)
% DATAFOLDAPPEND - please describe this function.
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

% Remember that black boundary lines are generated as bottom, left, top,
% right associated with i = 1,2,3,4, respectively
    
if alert == 0
    
    % Include special case for the first block
    if index == 1
        
        for i = 1:4

            % Assign null values for purpose of DXF Generation
            Struct1(i).x = [];
            Struct1(i).y = [];

        end
        
        for i = 5:size(Struct1, 2)
            
            % No mod term needed

            % Do not consider the m offset here, since it will always be 0
            Struct1(i).x = Struct1(i).x;
            Struct1(i).y = Struct1(i).y;

            % Plot each segment with offset factored in
            plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)
            
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

        infostruct(index).lmaxnet = infostruct(index).lmax + ...
            infostruct(index-1).lmaxnet;
        
        for i = 1:4

            % Assign null values for purpose of DXF Generation
            Struct1(i).x = [];
            Struct1(i).y = [];

        end

        % First fold
        for i = 5:size(Struct1, 2)

            % Factor in offset for every location along the structure.
            Struct1(i).x = Struct1(i).x + msum;
            Struct1(i).y = Struct1(i).y + infostruct(index-1).lmax + lmax_sum;

            % Plot each segment as before, but with offset factored in
            plot(a, Struct1(i).x, Struct1(i).y, 'color', Struct1(i).color)

        end

        lmax_sum = infostruct(index-1).lmax + lmax_sum;

    end    

end

end