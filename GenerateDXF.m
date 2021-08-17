% GenerateDXF
% Lasted edited ??/??/???? by ??

function GenerateDXF(filename, data)

% Create DXF file with data struct: 
% {x(colume),y(colume),cut(boolean),color([0,0,0])}
% and save as 'filename'

fid = dxf_open(filename);
for q = 1:length(data)
    % set color
    fid = dxf_set(fid,'Color',data(q).color);
    % plot polyline
    dxf_polyline(fid,data(q).x,data(q).y,zeros(size(data(q).x)));
end

dxf_close(fid);
end