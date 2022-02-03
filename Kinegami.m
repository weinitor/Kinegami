% Kinegami
% Last Edited 7/22/2021 by Lucien Peach

function [infostruct, TransformStruct, DataNet] = Kinegami(D, r, n, ...
    JointStruct, mirror, triple, theta_mod, fingertip, selfassign, ...
    TransformStruct, DXF, split, segmentation, plotoption)

    addpath('DXFLib_v0.9.1')

    % Initialize infostruct  
    % Determine whether waypoint integration will take place
    if strcmp(selfassign, 'true') == 1 
        num = 2 + 5*(size(JointStruct, 2) - 1);
        infostruct(num) = struct();
        init_size = num;
    else
        num = 2 + 15*(size(JointStruct, 2)-1);
        infostruct(num) = struct();
        init_size = num;
    end
    
    N = size(JointStruct, 2) - 1;
    
    % Tube
    tube_height = 0.1;
    
    [ls] = Origami_Tube_Parameters(r, n);

    infostruct(1).ls = ls;
    infostruct(1).r = r;

    % Outputs default tube parameters
    [dataFoldTube, m, lmax] = Origami_Tube_CreasePattern(n, ls, tube_height, r);

    infostruct(1).m = m;
    infostruct(1).lmax = lmax;
    infostruct(1).n = n;
    infostruct(1).type = dataFoldTube;
    infostruct(1).name = "Tube";
    
    if strcmp(selfassign, 'true') == 1
        
        % Run if Joint Assignment has been pre-assigned
        [TransformStruct] = SelfAssign(TransformStruct, r, n, JointStruct, N, plotoption); 
        
    else
        
        % Joint Placement and Planar Analysis
        [TransformStruct] = JointPlacement(D, r, n, JointStruct, N, fingertip, plotoption);
    
%         % Joint Assignment and Sphere Analysis for DH specs
%         [TransformStruct] = JointAssignment(D, r, n, JointStruct, N, theta_mod, fingertip, plotoption);
        
    end

    % Define h1 and h2
    for i = 1:init_size
        
        infostruct(i).h1 = 0;
        infostruct(i).h2 = 0;
        
    end
    
    % Check a vectors for optimal orientation
    for i = 1:N
        
        % Extract a and c vectors for frame being adjusted
        a = TransformStruct(i).Oc(:, 1);
        c = TransformStruct(i).Oc(:, 3);
        
        % Find distance vector between initial point and subsequent
        Vd = TransformStruct(i+1).Oc(:, 4) - TransformStruct(i).Oc(:, 4); 
        
        % If dot product is less than 0, flip direction of a
        if dot(Vd, a) < 0
            
            TransformStruct(i).Oc(:, 1) = -a;
            
            % Also flip c to keep consistent with RHR
            TransformStruct(i).Oc(:, 3) = -c;
            
        end
        
    end
    
    % Populate proximal and distal potential frames prior to looping
    for i = 1:N+1
    
        % Establish adjustment parameters
        TransformStruct(i).adjust = TransformStruct(i).rs * ...
            TransformStruct(i).Oc(:, 1);

        % Create new fields for editing
        TransformStruct(i).Op = TransformStruct(i).Oc;
        TransformStruct(i).Od = TransformStruct(i).Oc;

        % Edit new fields to create proximal and distal frames
        TransformStruct(i).Op(:, 4) = TransformStruct(i).Op(:, 4) ...
            - TransformStruct(i).adjust;
        TransformStruct(i).Od(:, 4) = TransformStruct(i).Od(:, 4) ...
            + TransformStruct(i).adjust;
        
    end
    
    % This figure call has no purpose except to prevent previous figure
    % from closing
    if strcmp(plotoption, 'on') == 1
        figure()
    end
    
    % Execute for JointPlacement Iteration
    if strcmp(selfassign, 'false') == 1
        
        % Create new structure for Joint & Waypoint Data Storage
        CumulativeStruct((3*N)+1) = struct();
        
        % Populating CumulativeStruct for main joints
        index = 0;
        for i = 1:3:(3*N)+1

            index = index+1;
            
            % Populate information for final joint plotted (initial joint
            % in series). No waypoint data here.
            CumulativeStruct(i).Op = TransformStruct(index).Op;
            CumulativeStruct(i).Od = TransformStruct(index).Od;
            
            % Converting JointStruct to CumulativeStruct
            CumulativeStruct(i).qm = JointStruct(index).qm;
            CumulativeStruct(i).q0 = JointStruct(index).q0;
            CumulativeStruct(i).type = JointStruct(index).type;
            CumulativeStruct(i).nz = JointStruct(index).nz;
            
            % Determine indexing for each of the main joints
            mainjoint = (index-1)*15+2;                
            
            % Classification by Joint Type
            if CumulativeStruct(i).type == 'R'

                theta_m = CumulativeStruct(i).qm;

                % Revolute Joint       
                [lengths, ls] = Origami_RevoluteJoint_Parameters(r, n, theta_m);

                infostruct(mainjoint).r = r;
                infostruct(mainjoint).ls = ls;
                infostruct(mainjoint).nz = CumulativeStruct(i).nz;
                nz = infostruct(mainjoint).nz;

                [dataFoldD, m, lmax] = Origami_RevoluteJoint_CreasePattern(lengths, ls, n, ...
                    infostruct(mainjoint).h1, infostruct(mainjoint).h2, r, theta_m, nz);

                infostruct(mainjoint).m = m;
                infostruct(mainjoint).lmax = lmax;
                infostruct(mainjoint).n = n;
                infostruct(mainjoint).type = dataFoldD; 
                infostruct(mainjoint).name = "Revolute";

            elseif CumulativeStruct(i).type == 'P'

                % Prismatic Joint
                beta = pi/3;
                nl = 2;
                d0 = CumulativeStruct(i).q0;

                [ls, l1, h0, dm, PJ_alpha] = Origami_PrismaticJoint_Parameters(r, n, beta, d0, nl);

                infostruct(mainjoint).r = r;
                infostruct(mainjoint).ls = ls;
                infostruct(mainjoint).n = n;
                infostruct(mainjoint).l1 = l1;
                infostruct(mainjoint).PJ_alpha = PJ_alpha;

                [dataFoldE, m, lmax] = Origami_PrismaticJoint_CreasePattern(r, n, nl, ls, l1, dm, h0, ...
                    infostruct(mainjoint).h1, infostruct(mainjoint).h2, PJ_alpha, beta);

                infostruct(mainjoint).m = m;
                infostruct(mainjoint).lmax = lmax;
                infostruct(mainjoint).type = dataFoldE;
                infostruct(mainjoint).name = "Prismatic";

            elseif CumulativeStruct(i).type == 'F'

                % Fingertip "Joint"
                theta_m = CumulativeStruct(i).qm;

                % Use same calculation file as revolute joint
                [lengths, ls] = Origami_RevoluteJoint_Parameters(r, n, theta_m);

                infostruct(mainjoint).r = r;
                infostruct(mainjoint).ls = ls;

                [FingertipFold, m, lmax] = Origami_Fingertip_CreasePattern(lengths, ls, n, ...
                    infostruct(mainjoint).h1, r, theta_m);

                infostruct(mainjoint).m = m;
                infostruct(mainjoint).lmax = lmax;
                infostruct(mainjoint).n = n;
                infostruct(mainjoint).type = FingertipFold; 
                infostruct(mainjoint).name = "Fingertip";
                
            end
            
        end
        
        index = 0;
        % Waypoint data
        for i = 2:3:(3*N)-1
            
            index = index+1;

            waypoint2index = (index-1)*15+7;
            waypoint1index = (index-1)*15+12;
            
            % Waypoint 2
            CumulativeStruct(i).type = 'W';
            CumulativeStruct(i).waypoint = TransformStruct(index+1).waypoint2;
            
            % Waypoint 2 Information 
            height = 0;
            
            [ls] = Origami_Tube_Parameters(r, n);
            infostruct(waypoint2index).ls = ls;
            infostruct(waypoint2index).r = r;
            
            % Outputs tube of height 0
            [dataFoldW, m, lmax] = Origami_Tube_CreasePattern(n, ls, height, r);
            
            infostruct(waypoint2index).m = m;
            infostruct(waypoint2index).lmax = lmax;
            infostruct(waypoint2index).n = n;
            infostruct(waypoint2index).type = dataFoldW;
            infostruct(waypoint2index).name = 'Waypoint';
            
            % Waypoint 1
            CumulativeStruct(i+1).type = 'W';
            CumulativeStruct(i+1).waypoint = TransformStruct(index+1).waypoint1;
            
            % Waypoint 1 Information
            height = 0;
            
            [ls] = Origami_Tube_Parameters(r, n);
            infostruct(waypoint1index).ls = ls;
            infostruct(waypoint1index).r = r;
            
            % Outputs tube of height 0
            [dataFoldW, m, lmax] = Origami_Tube_CreasePattern(n, ls, height, r);
            
            infostruct(waypoint1index).m = m;
            infostruct(waypoint1index).lmax = lmax;
            infostruct(waypoint1index).n = n;
            infostruct(waypoint1index).type = dataFoldW;
            infostruct(waypoint1index).name = 'Waypoint';
            
        end
        
        % At this point, we have populated the initial tube, all main
        % joints, and all waypoint joints. Now, we run DubinsTube to
        % connect these.
        index = 0;
        for i = 1:3:(3*N)-2
            
            % Indexing
            index = index+1;            
            val1 = (index-1)*15+3;
            val2 = (index-1)*15+8;
            val3 = (index-1)*15+13;

            % Run DubinsTube analysis
            [infostruct] = DubinsTube(r, n, CumulativeStruct(i).Od, ...
                CumulativeStruct(i+1).waypoint, infostruct, val1, mirror, split);
            [infostruct] = DubinsTube(r, n, CumulativeStruct(i+1).waypoint, ...
                CumulativeStruct(i+2).waypoint, infostruct, val2, mirror, split);
            [infostruct] = DubinsTube(r, n, CumulativeStruct(i+2).waypoint, ...
                CumulativeStruct(i+3).Op, infostruct, val3, mirror, split);
            
        end
        

    end
    
    % infostruct population for selfassign (no waypoints)
    if strcmp(selfassign, 'true') == 1
        for i = 1:N+1

            if JointStruct(i).type == 'R'

                theta_m = JointStruct(i).qm;
                jointindex = (i-1)*5+2;

                % Revolute Joint       
                [lengths, ls] = Origami_RevoluteJoint_Parameters(r, n, theta_m);

                infostruct(jointindex).r = r;
                infostruct(jointindex).ls = ls;
                infostruct(jointindex).nz = JointStruct(i).nz;
                nz = infostruct(jointindex).nz;

                [dataFoldD, m, lmax] = Origami_RevoluteJoint_CreasePattern(lengths, ls, n, ...
                    infostruct(i).h1, infostruct(i).h2, r, theta_m, nz);

                infostruct(jointindex).m = m;
                infostruct(jointindex).lmax = lmax;
                infostruct(jointindex).n = n;
                infostruct(jointindex).type = dataFoldD; 
                infostruct(jointindex).name = "Revolute";

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
                infostruct(jointindex).PJ_alpha = PJ_alpha;


                [dataFoldE, m, lmax] = Origami_PrismaticJoint_CreasePattern(r, n, nl, ls, l1, dm, h0, ...
                    infostruct(i).h1, infostruct(i).h2, PJ_alpha, beta);

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
                    infostruct(i).h1, r, theta_m);

                infostruct(jointindex).m = m;
                infostruct(jointindex).lmax = lmax;
                infostruct(jointindex).n = n;
                infostruct(jointindex).type = FingertipFold; 
                infostruct(jointindex).name = "Fingertip";
                
            else
                
                % Waypoint Specification
                jointindex = (i-1)*5+2;
                height = 0.001;

                [ls] = Origami_Tube_Parameters(r, n);

                infostruct(jointindex).ls = ls;
                infostruct(jointindex).r = r;

                % Outputs default tube parameters
                [dataFoldV, m, lmax] = Origami_Tube_CreasePattern(n, ls, height, r);

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

                [infostruct] = DubinsTube(r, n, TransformStruct(i).Od, ...
                    TransformStruct(i+1).Op, infostruct, newval, mirror, split);   
                
            end
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
            quiver3(TransformStruct(i).Od(1,4), TransformStruct(i).Od(2,4), ...
                TransformStruct(i).Od(3,4), TransformStruct(i).Od(1,j), ...
                TransformStruct(i).Od(2,j), TransformStruct(i).Od(3,j), ...
                'AutoScaleFactor', 0.025, 'Linewidth', 3, 'Color', colorvector(j));
            
            end
            
        end
        
        for j = 1:3

            % Plotting "distal" frame axes
            quiver3(TransformStruct(i).Op(1,4), TransformStruct(i).Op(2,4), ...
                TransformStruct(i).Op(3,4), TransformStruct(i).Op(1,j), ...
                TransformStruct(i).Op(2,j), TransformStruct(i).Op(3,j), ...
                'AutoScaleFactor', 0.025, 'Linewidth', 3, 'Color', colorvector(j));
        end
              
    end
    
    % Identify red, green, blue, as a, b, c
    legend('a', 'b', 'c')
    view(45, 45)
    daspect([1, 1, 1])
    
    % Add frames to previous plot
    for i = 1:N+1
        
        plot1 = frameplot(TransformStruct(i).Op, 'black'); % Plotting "proximal"
        
        if i < N+1 % again, do not worry about fingertip distal
            
            plot2 = frameplot(TransformStruct(i).Od, 'black'); % Plotting "distal"
        
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
        proximalcenter(i, 1) = TransformStruct(i).Od(1, 4);
        proximalcenter(i, 2) = TransformStruct(i).Od(2, 4);
        proximalcenter(i, 3) = TransformStruct(i).Od(3, 4);
        
        distalcenter(i, 1) = TransformStruct(i).Op(1, 4);
        distalcenter(i, 2) = TransformStruct(i).Op(2, 4);
        distalcenter(i, 3) = TransformStruct(i).Op(3, 4);
        
%         x = [proximalcenter(i, 1), distalcenter(i+1, 1)];
%         y = [proximalcenter(i, 2), distalcenter(i+1, 2)];
%         z = [proximalcenter(i, 3), distalcenter(i+1, 3)];
        
%         h = plot3(x, y, z, 'color', 'k', 'Linewidth', 1.5);
%         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
 
    end
    
    % WEI: Fixed the orientation of inital tube 1/27/2022
    % Add initial tube plotting. Take all 3 possible orientations into
    % account for accurate point plotting.
    
    % a = [1; 0; 0]
    if abs(TransformStruct(1).Op(1, 1)) == 1
        
        x = [distalcenter(1, 1), distalcenter(1, 1)]...
            -[sign(TransformStruct(1).Op(1, 1))*infostruct(1).lmax, 0];
        y = [distalcenter(1, 2), distalcenter(1, 2)];
        z = [distalcenter(1, 3), distalcenter(1, 3)];
        
        h = plot3(x, y, z, 'color', 'k', 'Linewidth', 4);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
    % a = [0; 1; 0]   
    elseif abs(TransformStruct(1).Op(2, 1)) == 1
        
        x = [distalcenter(1, 1), distalcenter(1, 1)];
        y = [distalcenter(1, 2), distalcenter(1, 2)]...
            -[sign(TransformStruct(1).Op(2, 1))*infostruct(1).lmax, 0];
        z = [distalcenter(1, 3), distalcenter(1, 3)];
        
        h = plot3(x, y, z, 'color', 'k', 'Linewidth', 4);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
    % a = [0; 0; 1]    
    else
        
        x = [distalcenter(1, 1), distalcenter(1, 1)];
        y = [distalcenter(1, 2), distalcenter(1, 2)];
        z = [distalcenter(1, 3), distalcenter(1, 3)]...
            -[sign(TransformStruct(1).Op(3, 1))*infostruct(1).lmax, 0];
        
        h = plot3(x, y, z, 'color', 'k', 'Linewidth', 4);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
    end
    
    % Colorvector data for joint spheres
    black = [0,0,0];
    red = [0.6350, 0.0780, 0.1840];
    orange = [0.8500, 0.3250, 0.0980];
    yellow = [0.929, 0.694, 0.125];
    green = [0.466, 0.674, 0.188];
    blue = [0, 0.447, 0.741];
    purple = [0.4940, 0.1840, 0.5560];

    colorvector = [black; red; orange; yellow; green; blue; purple];
    
    % Add joint spheres
    % For this specific instance of SphericalSampling, we require the
    % assignment of handle and thus must set plotoption to 'on' regardless
    % of its assignment.
    plotoption = 'on';
    
    for i = 1:N+1
        
        % Run SphericalSampling on each joint
        [TransformStruct(i).dubinsplot, handle] = SphericalSampling(TransformStruct(i).Oc(:, 4), ...
            TransformStruct(i).rs, colorvector(i, :), plotoption);   
        
        % Turn off legend for appearance
        handle.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
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
    
    % WEI: TO BE FIXED (take spliting into consideration)
    % Dubins Plotting    
    for i = 1:N
        
        val = 1 + (i-1)*5; 
        
        newval = val+2;
    
        [TransformStruct] = DubinsPlot(TransformStruct, infostruct, newval, i);
    
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
