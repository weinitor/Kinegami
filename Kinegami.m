% Kinegami
% Last Edited 6/28/2021 by Lucien Peach

function [infostruct, TransformStruct] = Kinegami(D, r, n, JointStruct, mirror, triple, theta_mod, fingertip)

    % Initialize infostruct
    num = 1 + 5*(size(JointStruct, 2) - 1);
    infostruct(num) = struct();
    
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
    
    % Joint Assignment and Sphere Analysis
    [TransformStruct] = JointAssignment(D, r, n, JointStruct, N, theta_mod, fingertip);

    % Define h1 and h2
    h1 = 0;
    h2 = 0;
    
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
    
    for i = 1:N
        
        if JointStruct(i).type == 'R'
            
            theta_m = JointStruct(i).qm;
            jointindex = (i-1)*5+2;
            
            % Revolute Joint
            [lengths, ls] = D_creasedesign_updated(r, n, theta_m);
            
            infostruct(jointindex).r = r;
            infostruct(jointindex).ls = ls;
            
            [dataFoldD, m, lmax] = D_papercut(lengths, ls, n, h1, h2, r, theta_m);
            
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
            
            
            [dataFoldE, m, lmax] = E_papercut(r, n, nl, ls, l1, dm, h0, h1, h2, E_alpha, beta);

            infostruct(jointindex).m = m;
            infostruct(jointindex).lmax = lmax;
            infostruct(jointindex).type = dataFoldE;
            infostruct(jointindex).name = "Prismatic";
            
        else
            
        end
        
        val = 1 + (i-1)*5; 
        
        newval = val+2;

        % Run Dubins Tube Analysis
        [infostruct] = DubinsTube(r, n, TransformStruct(i).Od, ...
            TransformStruct(i+1).Op, infostruct, newval, mirror);
      
    end
    
    % Add field for tracking lmaxnet
    infostruct(1).lmaxnet = infostruct(1).lmax;
    
    msum = 0;
    lmax_sum = 0;
    
    figure()
    set(gcf, 'color', 'w')
    hold on

    % Loop through indices to plot 
    for index = 1:size(infostruct, 2)

        [msum, lmax_sum, infostruct] = DataFoldAppend(infostruct(index).type, ...
            infostruct, index, msum, lmax_sum);

    end
    
    msum = 0;
    lmax_sum = 0;

    % Loop through indices to plot duplication
    for index = 1:size(infostruct, 2)

        [msum, lmax_sum] = DataFoldDuplicate(infostruct(index).type, ...
            infostruct, index, msum, lmax_sum, triple);

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
    
    % Run JointAssignment again for visualization
    [TransformStruct] = JointAssignment(D, r, n, JointStruct, N, theta_mod, fingertip);
    
end