function [infostruct, TransformStruct, DataNet, JointStruct] = Kinegami(D, r, n, ...
    JointStruct, mirror, triple, theta_mod, fingertip, ...
    TransformStruct, DXF, split, segmentation, plotoption, jointselect, ...
    tubeinit)
% KINEGAMI - Generates a crease pattern that folds into a serial mechanism
% from a given D-H specification.

% Inputs:
%   D               - the D-H parameter table of values: i x [a, alpha, d,
%                     theta].
%   r               - desired radius of folded origami linkage.
%   n               - number of sides of folded origami linkage.
%   JointStruct     - a data structure that contains information about
%                     joint parameters, frames, and connection pathways.
%   mirror          - setting required for creation of elbow schematic. See
%                     Origami_Elbow_CreasePattern.m for further detail.
%   triple          - a string input which controls the degree to which
%                     duplication takes place.
%   theta_mod       - revolute joint parameter for use within
%                     RotationalMatrix.m.
%   fingertip       - string input ('x', 'y', or 'z') used for fingertip
%                     orientation assignment.
%   TransformStruct - data structure which contains information about
%                     bounding spheres, associated planes, and other
%                     related information.
%   DXF             - string input which dicates DXF generation. 
%   split           - setting required for creation of elbow schematic. See
%                     Origami_Elbow_CreasePattern.m for further detail. 
%   segmentation    - string input which dictates segment splitting for
%                     ease of fabrication.
%   plotoption      - string input which dictates plotting.
%   jointselect     - string input which dictates type of joint placement.
%   tubeinit        - string input which dictates initial tube plotting

% Outputs:
%   infostruct      - updated data structure that includes all the
%                     information needed for the construction of the full
%                     schematic.
%   TransformStruct - updated data structure that contains information
%                     about bounding spheres, associated planes, and other
%                     related information.
%   DataNet         - cumulative data structure used for DXF generation.
%   JointStruct     - updated data structure that contains information
%                     about joint parameters, frames, and connection
%                     pathways.

% Authors: 
% Lucien Peach <peach@seas.upenn.edu>
% Last Edited 1/25/2023
%
% Copyright (C) 2022 The Trustees of the University of Pennsylvania. 
% All rights reserved. Please refer to LICENSE.md for detail.


    % Add all the folders and subfolders to the search path
    addpath(genpath(fileparts(mfilename('fullpath'))));
    
    % Define h1 and h2 value for extended revolute joints
    h1 = 0.035;
    h2 = 0.035;
    
    % Determine N (this will ultimately change for JointPlacement.m case)
    N = size(JointStruct, 2) - 1;
    
    % Define h1 and h2 data for each joint (in JointStruct)
    for i = 1:N+1
        
        if JointStruct(i).type == 'E'
            JointStruct(i).h1 = h1;
            JointStruct(i).h2 = h2;
        else
            JointStruct(i).h1 = 0;
            JointStruct(i).h2 = 0;
        end        
    end
    
    if strcmp(jointselect, 'selfassign') == 1
               
        % Run if Joint Assignment has been pre-assigned
        [TransformStruct] = SelfAssign(TransformStruct, r, n, JointStruct, N, plotoption); 
        
    else
        
        if strcmp(jointselect, 'placementB') == 1
            
            % Joint Placement and Planar Analysis
            [TransformStruct, JointStruct, N] = JointPlacementB(D, r, n, ...
                JointStruct, N, theta_mod, fingertip, plotoption);
            
        elseif strcmp(jointselect, 'placementA') == 1
    
            % Joint Assignment and Sphere Analysis for DH specs
            [TransformStruct] = JointPlacementA(D, r, n, JointStruct, N, ...
                theta_mod, fingertip, plotoption);
            
        end
        
    end
    
    % Initialize infostruct
    num = 2 + 5*(size(JointStruct, 2) - 1);
    infostruct(num) = struct();
    init_size = num;
    
    % Include initial tube to 0.1m in instance of tubeinit
    if strcmp(tubeinit, 'on') == 1
        tube_height = 0.1;
    else
        tube_height = 0;
    end

    % Populate initial tube in infostruct
    [ls] = Origami_Tube_Parameters(r, n);

    infostruct(1).ls = ls;
    infostruct(1).r = r;

    % Outputs default tube parameters
    [dataFoldTube, m, lmax] = Origami_Tube_CreasePattern(n, ls, tube_height);

    infostruct(1).m = m;
    infostruct(1).lmax = lmax;
    infostruct(1).n = n;
    infostruct(1).type = dataFoldTube;
    infostruct(1).name = "Tube";

    % Define h1 and h2
    for i = 1:init_size
        
        infostruct(i).h1 = 0;
        infostruct(i).h2 = 0;
        
    end
        
    % Populate proximal and distal potential frames prior to looping
    index = 0;
    for i = 1:N+1
        
        if strcmp(JointStruct(i).type, 'W') ~= 1
            
            index = index + 1;
    
            % Establish adjustment parameters
            TransformStruct(index).adjust = TransformStruct(index).rs * ...
                TransformStruct(index).Oc(:, 1);

            % Create new fields for editing
            JointStruct(i).Op = TransformStruct(index).Oc;
            JointStruct(i).Od = TransformStruct(index).Oc;

            % Edit new fields to create proximal and distal frames
            JointStruct(i).Op(:, 4) = JointStruct(i).Op(:, 4) ...
                - TransformStruct(index).adjust;
            JointStruct(i).Od(:, 4) = JointStruct(i).Od(:, 4) ...
                + TransformStruct(index).adjust;        
        end 
    end
    
    % This figure call has no purpose except to prevent previous figure
    % from closing
    if strcmp(plotoption, 'on') == 1
        figure()
    end
    
    % infostruct population and execution of DubinsTube()
    for i = 1:N+1

        if JointStruct(i).type == 'R' || JointStruct(i).type == 'E'
            
            jointindex = (i-1)*5+2;
            
            % for extended revolute, override h1 and h2 measurements
            % (hardcoded, define here)
            if JointStruct(i).type == 'E'
                % naming convention also included here
                infostruct(jointindex).name = "Extended Revolute";
                infostruct(jointindex).h1 = h1;
                infostruct(jointindex).h2 = h2;
            elseif JointStruct(i).type == 'R'
                infostruct(jointindex).name = "Revolute";                
            end

            theta_m = JointStruct(i).qm;

            % Revolute Joint       
            [lengths, ls] = Origami_RevoluteJoint_Parameters(r, n, theta_m);

            infostruct(jointindex).r = r;
            infostruct(jointindex).ls = ls;
            infostruct(jointindex).nz = JointStruct(i).nz;
            nz = infostruct(jointindex).nz;

            [dataFoldD, m, lmax] = Origami_RevoluteJoint_CreasePattern(lengths, ls, n, ...
                infostruct(jointindex).h1, infostruct(jointindex).h2, nz);

            infostruct(jointindex).m = m;
            infostruct(jointindex).lmax = lmax;
            infostruct(jointindex).n = n;
            infostruct(jointindex).type = dataFoldD; 

        elseif JointStruct(i).type == 'P'

            % Prismatic Joint
            beta = pi/3;
            nl = 2;
            d0 = JointStruct(i).q0;

            jointindex = (i-1)*5+2;

            [ls, l1, h0, dm, PJ_alpha] = Origami_PrismaticJoint_Parameters(r, n, beta, d0, nl);

            infostruct(jointindex).r = r;
            infostruct(jointindex).ls = ls;
            infostruct(jointindex).n = n;
            infostruct(jointindex).l1 = l1;
            infostruct(jointindex).h0 = h0;
            infostruct(jointindex).PJ_alpha = PJ_alpha;


            [dataFoldE, m, lmax] = Origami_PrismaticJoint_CreasePattern(n, nl, ls, l1, dm, ...
                infostruct(jointindex).h1, infostruct(jointindex).h2, PJ_alpha);

            infostruct(jointindex).m = m;
            infostruct(jointindex).lmax = lmax;
            infostruct(jointindex).type = dataFoldE;
            infostruct(jointindex).name = "Prismatic";

        elseif JointStruct(i).type == 'F'

            % Fingertip "Joint"
            theta_m = JointStruct(i).qm;
            jointindex = (i-1)*5+2;

            % Use same calculation file as revolute joint
            [lengths, ls] = Origami_RevoluteJoint_Parameters(r, n, theta_m);

            infostruct(jointindex).r = r;
            infostruct(jointindex).ls = ls;

            [FingertipFold, m, lmax] = Origami_Fingertip_CreasePattern(lengths, ls, n, ...
                infostruct(jointindex).h1);

            infostruct(jointindex).m = m;
            infostruct(jointindex).lmax = lmax;
            infostruct(jointindex).n = n;
            infostruct(jointindex).type = FingertipFold; 
            infostruct(jointindex).name = "Fingertip";

        else

            % Waypoint Specification
            jointindex = (i-1)*5+2;
            height = 0;

            [ls] = Origami_Tube_Parameters(r, n);

            infostruct(jointindex).ls = ls;
            infostruct(jointindex).r = r;

            % Outputs default tube parameters
            [dataFoldV, m, lmax] = Origami_Tube_CreasePattern(n, ls, height);

            infostruct(jointindex).m = m;
            infostruct(jointindex).lmax = lmax;
            infostruct(jointindex).n = n;
            infostruct(jointindex).type = dataFoldV;
            infostruct(jointindex).name = "Zero";
        end

        val = 1 + (i-1)*5; 

        newval = val+2;

        % Run Dubins Tube Analysis
        if i < N+1

            [infostruct] = DubinsTube(r, n, JointStruct(i).Od, ...
                JointStruct(i+1).Op, infostruct, newval, mirror, split);   
        end
    end
    
    % Create new figure to demonstrate the frames and their connections, as
    % well as the connection pipes
    figure()
    set(gcf, 'color', 'w')
    hold on
    
    % Indicate colors for a, b, c, for each frame
    colorvector = ['r'; 'g'; 'b'];
    
    for i = 1:N+1
        
        for j = 1:3
            
            % Do not worry about printing fingertip distal
            if i < N+1
            
            % Plotting "proximal" frame axes
            quiver3(JointStruct(i).Od(1,4), JointStruct(i).Od(2,4), ...
                JointStruct(i).Od(3,4), JointStruct(i).Od(1,j), ...
                JointStruct(i).Od(2,j), JointStruct(i).Od(3,j), ...
                'AutoScaleFactor', 0.025, 'Linewidth', 3, 'Color', colorvector(j));
            end     
        end
        
        for j = 1:3

            % Plotting "distal" frame axes
            quiver3(JointStruct(i).Op(1,4), JointStruct(i).Op(2,4), ...
                JointStruct(i).Op(3,4), JointStruct(i).Op(1,j), ...
                JointStruct(i).Op(2,j), JointStruct(i).Op(3,j), ...
                'AutoScaleFactor', 0.025, 'Linewidth', 3, 'Color', colorvector(j));
        end         
    end
    
    % Identify red, green, blue, as a, b, c
    legend('a', 'b', 'c')
    view(45, 45)
    daspect([1, 1, 1])
    
    % Add frames to previous plot
    for i = 1:N+1
        
        plot1 = frameplot(JointStruct(i).Op, 'black'); % Plotting "proximal"
        
        if i < N+1 % again, do not worry about fingertip distal
            
            plot2 = frameplot(JointStruct(i).Od, 'black'); % Plotting "distal"        
        end
        
        % Ignore in legend
        plot1.Annotation.LegendInformation.IconDisplayStyle = 'off';
        plot2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        % Recall that the first "proximal" is the distal of the first
        % joint. The first "distal" is the proximal of the second joint.      
    end
    
    % Initialize proximalcenter and distalcenter vectors (for plotting
    % frame connections)
    proximalcenter = zeros(N+1, 3);
    distalcenter = zeros(N+1, 3);
    
    % Add distal(i-1) - proximal(i) frame connections
    for i = 1:N
        
        % Identify location in space for "proximal" and "distal"
        proximalcenter(i, 1) = JointStruct(i).Od(1, 4);
        proximalcenter(i, 2) = JointStruct(i).Od(2, 4);
        proximalcenter(i, 3) = JointStruct(i).Od(3, 4);
        
        distalcenter(i, 1) = JointStruct(i).Op(1, 4);
        distalcenter(i, 2) = JointStruct(i).Op(2, 4);
        distalcenter(i, 3) = JointStruct(i).Op(3, 4);
        
%         x = [proximalcenter(i, 1), distalcenter(i+1, 1)];
%         y = [proximalcenter(i, 2), distalcenter(i+1, 2)];
%         z = [proximalcenter(i, 3), distalcenter(i+1, 3)];
        
%         h = plot3(x, y, z, 'color', 'k', 'Linewidth', 1.5);
%         h.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
    end
    
    % Add initial tube plotting. Take all 3 possible orientations into
    % account for accurate point plotting.
    
    if strcmp(tubeinit, 'on') == 1
        % a = [1; 0; 0]
        if abs(JointStruct(1).Op(1, 1)) == 1

            x = [distalcenter(1, 1), distalcenter(1, 1)]...
                -[sign(JointStruct(1).Op(1, 1))*infostruct(1).lmax, 0];
            y = [distalcenter(1, 2), distalcenter(1, 2)];
            z = [distalcenter(1, 3), distalcenter(1, 3)];

            h = plot3(x, y, z, 'color', 'k', 'Linewidth', 4);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';

        % a = [0; 1; 0]   
        elseif abs(JointStruct(1).Op(2, 1)) == 1

            x = [distalcenter(1, 1), distalcenter(1, 1)];
            y = [distalcenter(1, 2), distalcenter(1, 2)]...
                -[sign(JointStruct(1).Op(2, 1))*infostruct(1).lmax, 0];
            z = [distalcenter(1, 3), distalcenter(1, 3)];

            h = plot3(x, y, z, 'color', 'k', 'Linewidth', 4);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';

        % a = [0; 0; 1]    
        else

            x = [distalcenter(1, 1), distalcenter(1, 1)];
            y = [distalcenter(1, 2), distalcenter(1, 2)];
            z = [distalcenter(1, 3), distalcenter(1, 3)]...
                -[sign(JointStruct(1).Op(3, 1))*infostruct(1).lmax, 0];

            h = plot3(x, y, z, 'color', 'k', 'Linewidth', 4);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';

        end
    end
    
    % Colorvector data for joint spheres
    black = [0,0,0];
    red = [0.6350, 0.0780, 0.1840];
    orange = [0.8500, 0.3250, 0.0980];
    yellow = [0.929, 0.694, 0.125];
    green = [0.466, 0.674, 0.188];
    blue = [0, 0.447, 0.741];
    purple = [0.4940, 0.1840, 0.5560];

    colorvector = [black; red; orange; yellow; green; blue; purple; black];
    
    % Add joint spheres
    % For this specific instance of SphericalSampling, we require the
    % assignment of handle and thus must set plotoption to 'on' regardless
    % of its assignment.
    plotoption = 'on';
    
    index = 0;
    for i = 1:N+1
        
        if strcmp(JointStruct(i).type, 'W') ~= 1
            
            index = index + 1;
        
            % Run SphericalSampling on each joint
            [TransformStruct(index).dubinsplot, handle] = SphericalSampling(TransformStruct(index).Oc(:, 4), ...
                TransformStruct(index).rs, colorvector(index, :), plotoption);   

            % Turn off legend for appearance
            handle.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        end  
    end
    
    % Plot Settings
    grid on
    
    % Change as needed
%     xlim([-0.1, 0.2])
%     ylim([-0.1, 0.2])
%     zlim([-0.05, 0.25])
        
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Frame Connections')
    
    % Dubins Plotting    
    for i = 1:N
        
        val = 1 + (i-1)*5; 
        
        newval = val+2;
    
        [JointStruct] = DubinsPlot(JointStruct, infostruct, newval, i);  
    end
    
    % Add field for tracking lmaxnet
    infostruct(1).lmaxnet = infostruct(1).lmax;
    
    % Create duplicate structure fields for structure storage
    for i = 1:init_size
        
        infostruct(i).duplicate = infostruct(i).type;
        
        % Create triplet structure which will only be used in the instance
        % of a 'triplet' input
        infostruct(i).triple = infostruct(i).type;         
    end
    
    msum = 0;
    lmax_sum = 0;
    
    fig1 = figure();
    set(gcf, 'color', 'w')
    
    % 'renderer' ensures that  the file will have the ability to be
    % converted to an SVG upon completion
    set(fig1, 'renderer', 'painters');
    a = axes;
    hold on

    % Loop through indices to plot 
    for index = 1:size(infostruct, 2)

        [msum, lmax_sum, infostruct, infostruct(index).type] = DataFoldAppend(a, infostruct(index).type, ...
            infostruct, index, msum, lmax_sum);
    end
    
    msum = 0;
    lmax_sum = 0;

    % Loop through indices to plot duplication
    for index = 1:size(infostruct, 2)

        [msum, lmax_sum, infostruct, infostruct(index).duplicate, ...
            infostruct(index).triple] ...
            = DataFoldDuplicate(a, infostruct(index).duplicate, ...
            infostruct(index).triple, infostruct, index, msum, lmax_sum, triple);
    end
        
    for i = 1:size(infostruct, 2)
        
        if i == 1
            lmaxtotal = infostruct(i).lmax;
        else
            lmaxtotal = infostruct(i).lmax + lmaxtotal;
        end      
    end
    
%     totalsize = ((size(infostruct, 2) - 1) / 5);
%     
%     count = 1;
%     
    val = size(infostruct, 2);
    
    % Add boundary boxes around individual components, split by tubes
    for i = 1:size(infostruct, 2)
        
        if strcmp(infostruct(i).name, "Tube") == 1 && i ~= 1 ...
                || strcmp(infostruct(i).name, "Fingertip") == 1 && i ~= 1
            
            if strcmp(segmentation, "on") == 1

                [dataFoldDissect] = DissectPlot(n, ls, i, infostruct, val);

                infostruct(end+1).type = dataFoldDissect;
            
            end            
        end 
    end
    
    % Adding boundary box around entirety of printed region
    [dataFoldBoundary] = BoundaryPlot(n, ls, lmaxtotal);
    
    infostruct(end+1).type = dataFoldBoundary;
    
    axis off
    daspect([1, 1, 1])
    
    % Generate DXF Structure Creation and Generation
    structsize = size(infostruct, 2);
    
    for i = 1:init_size
        
        % Assign a field noted as size to identify structure scope
        infostruct(i).size = size(infostruct(i).type, 2);
        
        % sizecounter will be used to initialize DataNet
        if i == 1
            
            sizecounter = infostruct(i).size;
            
        else
            
            sizecounter = sizecounter + infostruct(i).size;            
        end   
    end
    
    % Assign .size field to additional entries
    for i = init_size:structsize
        
        infostruct(i).size = size(infostruct(i).type, 2);               
    end
    
    % Specify length of structure
    if strcmp(triple, 'triple') == 1
        mult = 3;
    else
        mult = 2;
    end
    
    % initialize for speed 
    DataNet(mult*sizecounter) = struct();
    offset = 0; % Net counter
    
    % Consolidating the information for GenerateDXF into a single
    % structure. 
    for i = 1:init_size
                
        for j = 1:infostruct(i).size
            
            DataNet(j+offset).x = infostruct(i).type(j).x;
            DataNet(j+offset).y = infostruct(i).type(j).y;
            DataNet(j+offset).color = infostruct(i).type(j).color;
            
            DataNet(j+offset+infostruct(i).size).x = infostruct(i).duplicate(j).x;
            DataNet(j+offset+infostruct(i).size).y = infostruct(i).duplicate(j).y;
            DataNet(j+offset+infostruct(i).size).color = infostruct(i).duplicate(j).color;
            
            % If three sections are desired
            if strcmp(triple, 'triple') == 1
                
                DataNet(j+offset+2*infostruct(i).size).x = infostruct(i).triple(j).x;
                DataNet(j+offset+2*infostruct(i).size).y = infostruct(i).triple(j).y;
                DataNet(j+offset+2*infostruct(i).size).color = infostruct(i).triple(j).color;               
            end              
        end
        
        offset = offset + mult*infostruct(i).size;
       
    end
    
    % Adding the outline boxes, which only need to be printed once
    for i = init_size+1:structsize
        
        for j = 1:infostruct(i).size
            
            DataNet(j+offset).x = infostruct(i).type(j).x;
            DataNet(j+offset).y = infostruct(i).type(j).y;
            DataNet(j+offset).color = infostruct(i).type(j).color;        
        end      
        
        offset = offset + infostruct(i).size;     
    end
    
    % Set axis limits so that only print area is displayed
%     axis([0, (n+1)*ls, 0, infostruct(val).lmaxnet])
    
    % Implement feature so that data outside of print area is covered
    fillx1 = [-2*(n+1)*ls, 0, 0, -2*(n+1)*ls];
    filly1 = [0, 0, infostruct(val).lmaxnet + 0.01, infostruct(val).lmaxnet + 0.01];
    patch(fillx1, filly1, 'w', 'EdgeColor', 'none');
    
    fillx2 = [(n+1)*ls, 4*(n+1)*ls, 4*(n+1)*ls, (n+1)*ls];
    filly2 = [0, 0, infostruct(val).lmaxnet + 0.01, infostruct(val).lmaxnet + 0.01];
    patch(fillx2, filly2, 'w', 'EdgeColor', 'none');
    
    % If desired, generate a new DXF file
    if strcmp(DXF, 'on') == 1
        filename = 'KinegamiTest.dxf';
        GenerateDXF(filename, DataNet)
    end
    
end
