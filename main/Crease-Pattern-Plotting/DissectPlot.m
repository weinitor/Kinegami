function [DissectPlot] = DissectPlot(n, ls, index, infostruct, val)
% DISSECTPLOT - Creates boundary boxes around individual Dubins Linkage
% Segments, separated by interposed tube segments. Boundary range can be
% adjusted in Kinegami.m.

% Inputs:
%   n           - number of sides of folded origami linkage.
%   ls          - side length of folded origami linkage. 
%   index       - counting term which is incremented outside of this loop.
%                 For use in indexing of infostruct terms.
%   infostruct  - a data structure that includes all the information needed
%                 for the construction of the full schematic.
%   val         - value representative of the max size of the infostruct 
%                 data structure 

% Outputs:
%   DissectPlot - a data structure which contains pertinent information for
%                 the graphing of boundary boxes.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last Edited 8/8/2021
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


% Counter used for data structure indexing
count = 1;

% Identify color
black = [0, 0, 0];
blue = [0, 0, 1];

% Boundary coordinates

% Special case for initial scenario
if index == 5
    
    % Store initial tube information to boundary matrix
    boundary = [0, 0; ...
        (n+1)*ls, 0; ...
        (n+1)*ls, infostruct(index-1).lmaxnet + infostruct(index).lmax/2; ...
        0, infostruct(index-1).lmaxnet + infostruct(index).lmax/2; ...
        0, 0];
    
    % Store to data structure and plot
    DissectPlot(count).x = boundary(:, 1);
    DissectPlot(count).y = boundary(:, 2);
    DissectPlot(count).color = black;
    
    plot(DissectPlot(count).x, DissectPlot(count).y, 'color', ...
        DissectPlot(count).color)
    
    % Increase counter
    count = count + 1;
    
    % Add overlap indicator line
    overlap = [0, infostruct(index-1).lmaxnet + infostruct(index).lmax/2; ...
        (n+1)*ls, infostruct(index-1).lmaxnet + infostruct(index).lmax/2];
    
    % Store to structure and plot
    DissectPlot(count).x = overlap(:, 1);
    DissectPlot(count).y = overlap(:, 2);
    DissectPlot(count).color = blue;
    
    plot(DissectPlot(count).x, DissectPlot(count).y, 'color', ...
        DissectPlot(count).color)
    
    % Plot settings
    set(gcf, 'color', 'w')
    daspect([1, 1, 1])
    axis off
    
elseif index ~=5 && index ~= val
    
    % Address top and bottom height parameters for easy translation into
    % the boundary matrix
    bottomheight = infostruct(index-6).lmaxnet + infostruct(index-5).lmax/2 ...
        - 0.02;
    topheight = infostruct(index-1).lmaxnet + infostruct(index).lmax/2;
    
    % Store subsequent tube information and boundary matrix
    boundary = [0, bottomheight; ...
        (n+1)*ls, bottomheight; ...
        (n+1)*ls, topheight; ...
        0, topheight; ...
        0, bottomheight];
    
    % Store to data structure and plot
    DissectPlot(count).x = boundary(:, 1);
    DissectPlot(count).y = boundary(:, 2);
    DissectPlot(count).color = black;
    
    plot(DissectPlot(count).x, DissectPlot(count).y, 'color', ...
        DissectPlot(count).color)
    
    % Increase counter
    count = count + 1;
    
    % Add overlap indicator line
    overlap = [0, topheight; ...
        (n+1)*ls, topheight];
    
    % Store to structure and plot
    DissectPlot(count).x = overlap(:, 1);
    DissectPlot(count).y = overlap(:, 2);
    DissectPlot(count).color = blue;
    
    plot(DissectPlot(count).x, DissectPlot(count).y, 'color', ...
        DissectPlot(count).color)
    
    % Plot settings
    set(gcf, 'color', 'w')
    daspect([1, 1, 1])
    axis off
    
else
        
    % Address top and bottom height parameters for easy translation into
    % the boundary matrix
    bottomheight = infostruct(index-3).lmaxnet + infostruct(index-2).lmax/2 ...
        - 0.02;
    topheight = infostruct(index).lmaxnet;
    
    % Store subsequent tube information and boundary matrix
    boundary = [0, bottomheight; ...
        (n+1)*ls, bottomheight; ...
        (n+1)*ls, topheight; ...
        0, topheight; ...
        0, bottomheight];
    
    % Store to data structure and plot
    DissectPlot(count).x = boundary(:, 1);
    DissectPlot(count).y = boundary(:, 2);
    DissectPlot(count).color = black;
    
    plot(DissectPlot(count).x, DissectPlot(count).y, 'color', ...
        DissectPlot(count).color)
    
    % Plot settings
    set(gcf, 'color', 'w')
    daspect([1, 1, 1])
    axis off
    hold on
    
end 

end