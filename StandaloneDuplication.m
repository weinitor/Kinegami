function [DataFoldNew] = StandaloneDuplication(DataFoldOld, ls, n, lmax, type, h1, maxlengths)
% STANDALONEDUPLICATION - The aim of this program will be to add the
% additional frame to each xxx_xxx_Standalone.m file, which is needed for
% overlap.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Created 12/13/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Add to right side for all types except 'revolute' (for ease of folding)
if strcmp(type, 'revolute') ~= 1
    
    % Determine total size of structure
    maxsize = size(DataFoldOld, 2);

    DataFoldNew(2*maxsize) = struct();

    % Create DataFoldNew
    for i = 1:maxsize

        DataFoldNew(i).x = DataFoldOld(i).x;
        DataFoldNew(i).y = DataFoldOld(i).y;
        DataFoldNew(i).color = DataFoldOld(i).color;

    end

    % Create Duplicate Entries
    for j = (maxsize+1):(2*maxsize)

        z = j - maxsize;

        DataFoldNew(j).x = DataFoldOld(z).x + n*ls;
        DataFoldNew(j).y = DataFoldOld(z).y;
        DataFoldNew(j).color = DataFoldOld(z).color;

    end

    % Erase Double Lines
    for k = (maxsize+1):(maxsize+4)

        DataFoldNew(k).x = [];
        DataFoldNew(k).y = [];
        DataFoldNew(k).color = [];

    end

    for k = 1:4

        DataFoldNew(k).x = [];
        DataFoldNew(k).y = [];
        DataFoldNew(k).color = [];

    end

    % Replace initial 4 entries (old boundary) with new boundary
    DataFoldNew(end+1).x = [0; (n+1)*ls];
    DataFoldNew(end).y = [0; 0];
    DataFoldNew(end+1).x = [(n+1)*ls; (n+1)*ls];
    DataFoldNew(end).y = [0; lmax];
    DataFoldNew(end+1).x = [(n+1)*ls; 0];
    DataFoldNew(end).y = [lmax; lmax];
    DataFoldNew(end+1).x = [0; 0];
    DataFoldNew(end).y = [lmax; 0];

    % Note that sequential (end+1) ---> (end) indexing due to new (end) value

    for i = 3:-1:0

        DataFoldNew(end-i).color = [0, 0, 0];

    end

    % Eliminate Printing Past Boundaries
    for index = 1:(2*maxsize)

        arraysizex = size(DataFoldNew(index).x, 1);

        for num = 1:arraysizex

            if DataFoldNew(index).x(num) > ((n+1)*ls)

                DataFoldNew(index).x(num) = ((n+1)*ls);

            end

        end

    end

end

% Determine duplication pattern for type: 'revolute'
if strcmp(type, 'revolute') == 1
    
    % Add duplication frame to LEFT of initial framing
    maxsize = size(DataFoldOld, 2);

    DataFoldNew(maxsize) = struct();

    % Create DataFoldNew
    for i = 1:maxsize

        DataFoldNew(i).x = DataFoldOld(i).x;
        DataFoldNew(i).y = DataFoldOld(i).y;
        DataFoldNew(i).color = DataFoldOld(i).color;

    end
    
    % Retool boundary limits
    DataFoldNew(1).x = [-ls; n*ls];
    DataFoldNew(1).y = [0; 0];
    DataFoldNew(2).x = [n*ls; n*ls];
    DataFoldNew(2).y = [0; lmax];
    DataFoldNew(3).x = [n*ls; -ls];
    DataFoldNew(3).y = [lmax; lmax];
    DataFoldNew(4).x = [-ls; -ls];
    DataFoldNew(4).y = [lmax; 0];
    
    % Add lines for consistent printing
    
    % Bottom Blue
    DataFoldNew(end+1).x = [-ls; 0];
    DataFoldNew(end).y = [h1; h1];
    DataFoldNew(end).color = [0, 0, 1];
    
    % Mid Orange
    DataFoldNew(end+1).x = [-ls; 0];
    DataFoldNew(end).y = [h1 + maxlengths; h1+maxlengths];
    DataFoldNew(end).color = [1, 0.41, 0]; 
    
    % Top Blue
    DataFoldNew(end+1).x = [-ls; 0];
    DataFoldNew(end).y = [h1 + 2*maxlengths; h1 + 2*maxlengths];
    DataFoldNew(end).color = [0, 0, 1];
    
end
    
% For instance of twist, which does not have initial vertical line, insert
% this vertical line into dataFoldNew
if strcmp(type, 'twist') == 1
    
    DataFoldNew(end+1).x = [n*ls, n*ls];
    DataFoldNew(end).y = [0; h1];
    DataFoldNew(end).color = [0, 0, 1];
    
end
    
end 