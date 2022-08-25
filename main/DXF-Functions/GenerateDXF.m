function GenerateDXF(filename, data)
% GENERATEDXF - Create DXF file with data struct:
% {x(colume),y(colume),cut(boolean),color([0,0,0])} and save as 'filename'.
% Uses the submodule "dxflib".

% Inputs:
%   filename    - name under which the created .dxf file will be saved
%   data        - data structure which will be used to the create the .dxf
%                 file. 

% Outputs:
%   None

% Authors:
% Unknown

fid = dxf_open(filename);
for q = 1:length(data)
    % set color
    fid = dxf_set(fid,'Color',data(q).color);
    % plot polyline
    dxf_polyline(fid,data(q).x,data(q).y,zeros(size(data(q).x)));
end

dxf_close(fid);
end