function [F, NDU, dzero,clampedNodes, contactNodes] = ApplyLoadsBoundary(Ft, Fc, Tc, botedge, rightcham, topredge, n_nodes, X, Y, Xb, Yb, X_body, Y_body, X_tri, Y_tri, bevelled, ANSYSmesh, normalMesh)
    F = zeros(2*n_nodes,1);

  
    %Find contact the nodes on the top surface between x = 0 and L_ref (user defined refine area)

    tolY = 1e-10;
    contactNodes = find( ...
        X >= 0.0000 & X <= Tc & abs(Y - 0.0000) < tolY);

    if isempty(contactNodes)
        error('No contact nodes found for contact line load.');
    end

  
    %Sort the contact nodes from left to right
    
    [~,idx] = sort(X(contactNodes));
    contactNodes = contactNodes(idx);

    disp('Contact nodes used for line load = ');
    disp(contactNodes);

    disp('Number of contact nodes = ');
    disp(length(contactNodes));

  
    %Solving for total contact length
    xContact = X(contactNodes);

    Lc = xContact(end) - xContact(1);

    if Lc <= 0
        error('Contact length is zero or negative.');
    end

   
    %Force per unit length solve
   
    qx = Fc / Lc;
    qy = Ft / Lc;


    %Tributary length for each node
    tributary = zeros(length(contactNodes),1);

    for i = 1:length(contactNodes)

        if i == 1
            tributary(i) = (xContact(i+1) - xContact(i)) / 2;

        elseif i == length(contactNodes)
            tributary(i) = (xContact(i) - xContact(i-1)) / 2;

        else
            tributary(i) = (xContact(i+1) - xContact(i-1)) / 2;
        end

    end


    %Convert the line loads into nodal forces
    Fx_nodes = qx * tributary;
    Fy_nodes = qy * tributary;

    %Check the totals
    disp('Total Fx applied = ');
    disp(sum(Fx_nodes));

    disp('Total Fy applied = ');
    disp(sum(Fy_nodes));

    %Build global load vector so assign proper x and y values
    for i = 1:length(contactNodes)

        n = contactNodes(i);

        F(2*n - 1) = F(2*n - 1) + Fx_nodes(i);
        F(2*n)     = F(2*n)     + Fy_nodes(i);

    end

%=================APPLICATION OF TOOL CLAMPING===========
    toll = 1e-6;
    dzero = [];


%===============FOR ANSYS MESH============================
    if ANSYSmesh == true
  % =====================================================
% BOTTOM EDGE (actually MIN Y)
% =====================================================
if botedge == true

    bottomNodes = find(abs(Y - min(Y)) < toll);

    dzero = [dzero;
             reshape([2*bottomNodes-1; 2*bottomNodes], [], 1)];

end


% =====================================================
% RIGHT CHAMFER EDGE (correct orientation)
% =====================================================
if rightcham == true

    % -------------------------------------------------
    % pick chamfer endpoints using TRUE extremes
    % (right side = largest X values)
    % -------------------------------------------------
    [~, i1] = max(X - Y);   % lower right-ish chamfer end
    [~, i2] = max(X + Y);   % upper right-ish chamfer end

    x1 = X(i1);  y1 = Y(i1);
    x2 = X(i2);  y2 = Y(i2);

    % line equation
    A = y2 - y1;
    B = x1 - x2;
    C = x2*y1 - x1*y2;

    dist = abs(A*X + B*Y + C) / sqrt(A^2 + B^2);

    den = (x2 - x1)^2 + (y2 - y1)^2;

    t = ((X - x1).*(x2 - x1) + (Y - y1).*(y2 - y1)) ./ den;

    rightEdgeNodes = find(dist < toll & t >= 0 & t <= 1);

    dzero = [dzero;
             reshape([2*rightEdgeNodes-1; 2*rightEdgeNodes], [], 1)];

end
    
% =========================================================
% TOP-RIGHT EDGE NODES (ANSYS mesh based)
% =========================================================
if topredge == true

    yMax = max(Y);
    xMin = min(X);
    xMax = max(X);

    xMid = (xMin + xMax) / 2;   % split domain into left/right

    topRightNodes = find( ...
        abs(Y - yMax) < toll & ...   % top boundary
        X >= xMid);                  % right half only

    for k = 1:length(topRightNodes)
        n = topRightNodes(k);
        dzero = [dzero;
                 2*n-1;
                 2*n];
    end

end
    end





 %If the tool is bevelled it runs this
    if bevelled == true

    if botedge == true
    bottomNodes = find(abs(Y - Y_body(1)) < toll);

        for k = 1:length(bottomNodes)
            n = bottomNodes(k);
            dzero = [dzero; 2*n-1; 2*n];
   
        end

     end


    if rightcham == true
    %Line equation coefficients
    x1r = X_body(3); 
    x2r = X_body(2); 
    y1r = Y_body(3);
    y2r = Y_body(2);
    A = y2r - y1r;
    B = x1r - x2r;
    C = x2r*y1r - x1r*y2r;
    
    dist = abs(A*X + B*Y + C) / sqrt(A^2 + B^2);
    
    %rightEdgeNodes = find(dist < toll);
    
    rightEdgeNodes = find(dist < toll);
        for k = 1:length(rightEdgeNodes)
            n = rightEdgeNodes(k);
            dzero = [dzero; 2*n-1; 2*n];
        end
    
    end
    
    if topredge == true
    %Top right half of top edge
        x_top_mid = (X_body(4) + X_body(3)) / 2;
        topRightNodes = find( ...
            X >= x_top_mid & ...
            X <= X_body(3) & ...
            abs(Y - Y_body(3)) < toll);
        for k = 1:length(topRightNodes)
            n = topRightNodes(k);
            dzero = [dzero; 2*n-1; 2*n];
        end
    
    end
    end





   if normalMesh == true


    % Bottom edge
    if botedge == true
    bottomNodes = find(abs(Y - Yb(1)) < toll);

        for k = 1:length(bottomNodes)
            n = bottomNodes(k);
            dzero = [dzero; 2*n-1; 2*n];
   
        end

    end


    if rightcham == true
    %Line equation coefficients
    x1r = Xb(3); 
    x2r = Xb(2); 
    y1r = Yb(3);
    y2r = Yb(2);
    A = y2r - y1r;
    B = x1r - x2r;
    C = x2r*y1r - x1r*y2r;
    
    dist = abs(A*X + B*Y + C) / sqrt(A^2 + B^2);
    
    %rightEdgeNodes = find(dist < toll);
    
    rightEdgeNodes = find(dist < toll);
        for k = 1:length(rightEdgeNodes)
            n = rightEdgeNodes(k);
            dzero = [dzero; 2*n-1; 2*n];
        end
    
    end
    
    if topredge == true
    %Top right half of top edge
        x_top_mid = (Xb(4) + Xb(3)) / 2;
        topRightNodes = find( ...
            X >= x_top_mid & ...
            X <= Xb(3) & ...
            abs(Y - Yb(3)) < toll);
        for k = 1:length(topRightNodes)
            n = topRightNodes(k);
            dzero = [dzero; 2*n-1; 2*n];
        end
    
    end
    end




    dzero = unique(dzero);
    NDU = length(dzero);
    clampedNodes = unique(ceil(dzero/2));

end