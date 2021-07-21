% Kinegami
% Last Edited 7/20/2021 by Lucien Peach

function [infostruct, TransformStruct, DataNet] = Kinegami(D, r, n, ...
    JointStruct, mirror, triple, theta_mod, fingertip, selfassign, ...
    TransformStruct)

    % Initialize infostruct
    num = 1 + 5*(size(JointStruct, 2) - 1);
    infostruct(num) = struct();
    init_size = num;
    
    N = size(JointStruct, 2) - 1;
    
    % Tube
    tube_height = 0.1;
    
    [ls] = Default_creasedesign(r, n);

    infostruct(1).ls = ls;
    infostruct(1).r = r;

    % Outputs default tube parameters
    [dataFoldTube, m, lmax] = Default_papercut(n, ls, tube_height, r);

    infostruct(1).m = m;
    infostruct(1).lmax = lmax;
    infostruct(1).n = n;
    infostruct(1).type = dataFoldTube;
    infostruct(1).name = "Tube";
    
    if strcmp(selfassign, 'true') == 1
        
        % Run if Joint Assignment has been pre-assigned
        [TransformStruct] = SelfAssign(TransformStruct, r, n, JointStruct, N); 
        
    else
    
        % Joint Assignment and Sphere Analysis for DH specs
        [TransformStruct] = JointAssignment(D, r, n, JointStruct, N, theta_mod, fingertip);
        
    end

    % Define h1 and h2
    for i = 1:init_size
        
        infostruct(i).h1 = 0;
        infostruct(i).h2 = 0;
        
    end
    
    % Check a vectors for optimal orientation
    for i = 1:N
        
        % Extract a and c vectors for frame being adjusted
        a = TransformStruct(i+1).Oc(:, 1);
        c = TransformStruct(i+1).Oc(:, 3);
        
        % Find distance vector between initial point and subsequent
        Vd = TransformStruct(i+1).Oc(:, 4) - TransformStruct(i).Oc(:, 4); 
        
        % If dot product is less than 0, flip direction of a
        if dot(Vd, a) < 0
            
            TransformStruct(i+1).Oc(:, 1) = -a;
            
            % Also flip c to keep consistent with RHR
            TransformStruct(i+1).Oc(:, 3) = -c;
            
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
        
    % Create new figure to demonstrate the frames and their connections, as
    % well as the connection pipes
    figure()
    set(gcf, 'color', 'w')
    hold on
    
    % Indicate colors for a, b, c, for each frame
    colorvector = ['r'; 'g'; 'b'];
    
    for i = 1:N
        
        for j = 1:3
            
            % Plotting "proximal" frame axes
            quiver3(TransformStruct(i).Od(1,4), TransformStruct(i).Od(2,4), ...
                TransformStruct(i).Od(3,4), TransformStruct(i).Od(1,j), ...
                TransformStruct(i).Od(2,j), TransformStruct(i).Od(3,j), ...
                'AutoScaleFactor', 0.025, 'Linewidth', 3, 'Color', colorvector(j));
            
        end
        
        for j = 1:3

            % Plotting "distal" frame axes
            quiver3(TransformStruct(i+1).Op(1,4), TransformStruct(i+1).Op(2,4), ...
                TransformStruct(i+1).Op(3,4), TransformStruct(i+1).Op(1,j), ...
                TransformStruct(i+1).Op(2,j), TransformStruct(i+1).Op(3,j), ...
                'AutoScaleFactor', 0.025, 'Linewidth', 3, 'Color', colorvector(j));
        end
              
    end
    
    % Identify red, green, blue, as a, b, c
    legend('a', 'b', 'c')
    daspect([1, 1, 1])
    
    % Add frames to previous plot
    for i = 1:N
        
        plot1 = frameplot(TransformStruct(i).Od, 'black'); % Plotting "proximal"
        plot2 = frameplot(TransformStruct(i+1).Op, 'black'); % Plotting "distal"
        
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
        
        distalcenter(i+1, 1) = TransformStruct(i+1).Op(1, 4);
        distalcenter(i+1, 2) = TransformStruct(i+1).Op(2, 4);
        distalcenter(i+1, 3) = TransformStruct(i+1).Op(3, 4);
        
        x = [proximalcenter(i, 1), distalcenter(i+1, 1)];
        y = [proximalcenter(i, 2), distalcenter(i+1, 2)];
        z = [proximalcenter(i, 3), distalcenter(i+1, 3)];
        
        h = plot3(x, y, z, 'color', 'k', 'Linewidth', 1.5);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
 
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
    
    % This figure call has no purpose except to prevent previous figure
    % from closing
    figure()
    for i = 1:N
        
        if JointStruct(i).type == 'R'
            
            theta_m = JointStruct(i).qm;
            jointindex = (i-1)*5+2;
            
            % Revolute Joint
            [lengths, ls] = D_creasedesign_updated(r, n, theta_m);
            
            infostruct(jointindex).r = r;
            infostruct(jointindex).ls = ls;
            
            [dataFoldD, m, lmax] = D_papercut(lengths, ls, n, ...
                infostruct(i).h1, infostruct(i).h2, r, theta_m);
            
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
            
            [ls, l1, h0, dm, E_alpha] = E_creasedesign(r, n, beta, d0, nl);
            
            infostruct(jointindex).r = r;
            infostruct(jointindex).ls = ls;
            infostruct(jointindex).n = n;
            infostruct(jointindex).l1 = l1;
            infostruct(jointindex).E_alpha = E_alpha;
            
            
            [dataFoldE, m, lmax] = E_papercut(r, n, nl, ls, l1, dm, h0, ...
                infostruct(i).h1, infostruct(i).h2, E_alpha, beta);

            infostruct(jointindex).m = m;
            infostruct(jointindex).lmax = lmax;
            infostruct(jointindex).type = dataFoldE;
            infostruct(jointindex).name = "Prismatic";
            
        else
            
        end
        
        val = 1 + (i-1)*5; 
        
        newval = val+2;

        % Run Dubins Tube Analysis
        [infostruct, tvec] = DubinsTube(r, n, TransformStruct(i).Od, ...
            TransformStruct(i+1).Op, infostruct, newval, mirror);     
      
    end
    
    % Dubins Plotting
    figure()
    set(gcf, 'color', 'w')
    hold on
    
    for i = 1:N
        
        val = 1 + (i-1)*5; 
        
        newval = val+2;
    
        [TransformStruct] = DubinsPlot(TransformStruct, infostruct, newval, i, tvec);
    
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
    
    figure()
    set(gcf, 'color', 'w')
    hold on

    % Loop through indices to plot 
    for index = 1:size(infostruct, 2)

        [msum, lmax_sum, infostruct, infostruct(index).type] = DataFoldAppend(infostruct(index).type, ...
            infostruct, index, msum, lmax_sum);

    end
    
    msum = 0;
    lmax_sum = 0;

    % Loop through indices to plot duplication
    for index = 1:size(infostruct, 2)

        [msum, lmax_sum, infostruct, infostruct(index).duplicate, ...
            infostruct(index).triple] ...
            = DataFoldDuplicate(infostruct(index).duplicate, ...
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
     val = size(infostruct, 2) - 1;
    
    
    for i = 1:size(infostruct, 2)
        
        if strcmp(infostruct(i).name, "Tube") == 1 && i ~= 1
            
%             count = count + 1;
            
            [dataFoldDissect] = DissectPlot(n, ls, i, infostruct, val);
            
            infostruct(end+1).type = dataFoldDissect;
            
        end
        
    end
    
    %    
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
    for i = init_size+1:structsize
        
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
    
    filename = 'KinegamiTest.dxf';
    GenerateDXF(filename, DataNet)
    
end