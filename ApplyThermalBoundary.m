function [fixedNodes, fixedValues, thermalTipNodes,ConvNodes] = ApplyThermalBoundary( ...
    thermalContact, ...
    botedge, rightcham, topredge, ...
    n_nodes, X, Y, ...
    Xb, Yb, ...
    X_body, Y_body, ...
    X_tri, Y_tri, ...
    bevelled, ANSYSmesh, normalMesh, ...
    Tambient, contactTemp)


% APPLY THERMAL BOUNDARY CONDITIONS (DIRICHLET ONLY)
% fixedNodes is nodes with prescribed temperature
% fixedValues istemperature values
fixedNodes  = [];
fixedValues = [];

tol = 1e-6;


%CONTACT REGION HOT TOOL INTERFACE
if bevelled == true
x1 = 0.0002533;  y1 = 0.0;
x2 = 0.0000123;  y2 = -0.0002121;

% Edge vector and length
dx = x2 - x1;
dy = y2 - y1;
L  = sqrt(dx^2 + dy^2);

%For each node compute 
%t = projection parameter along the edge
%d = perpendicular distance from the line
t = ((X - x1)*dx + (Y - y1)*dy) / L^2;
d = abs((X - x1)*dy - (Y - y1)*dx) / L;

tolD = 1e-9;%perpendicular distance tolerance
thermalTipNodes = find(d < tolD & t >= -tolD & t <= 1 + tolD);
else    
thermalTipNodes = find( ...
    X >= 0.0 & ...
    X <= thermalContact & ...
    abs(Y - 0.0) < tol);
end
if ~isempty(thermalTipNodes)

    %Apply contact temperature
    fixedNodes  = thermalTipNodes(:);
    fixedValues = contactTemp * ones(length(thermalTipNodes),1);

end


%ANSYS MESH BOUNDARIES
if ANSYSmesh == true

    if botedge == true

        bottomNodes = find(abs(Y - min(Y)) < tol);

        fixedNodes  = [fixedNodes; bottomNodes(:)];
        fixedValues = [fixedValues; Tambient * ones(length(bottomNodes),1)];

    end

    if rightcham == true

        [~, i1] = max(X - Y);
        [~, i2] = max(X + Y);

        x1 = X(i1); y1 = Y(i1);
        x2 = X(i2); y2 = Y(i2);

        A = y2 - y1;
        B = x1 - x2;
        C = x2*y1 - x1*y2;

        dist = abs(A*X + B*Y + C) / sqrt(A^2 + B^2);

        s = ((X - x1).*(x2 - x1) + ...
             (Y - y1).*(y2 - y1)) ./ ((x2 - x1)^2 + (y2 - y1)^2);

        rightNodes = find(dist < tol & s >= 0 & s <= 1);

        fixedNodes  = [fixedNodes; rightNodes(:)];
        fixedValues = [fixedValues; Tambient * ones(length(rightNodes),1)];

    end

    if topredge == true

        yTop = max(Y);
        xMid = (min(X) + max(X)) / 2;

        topNodes = find(abs(Y - yTop) < tol & X >= xMid);

        fixedNodes  = [fixedNodes; topNodes(:)];
        fixedValues = [fixedValues; Tambient * ones(length(topNodes),1)];

    end


%Bevelled tool geometry

elseif bevelled == true

    if botedge == true

        bottomNodes = find(abs(Y - min(Y)) < tol);

        fixedNodes  = [fixedNodes; bottomNodes(:)];
        fixedValues = [fixedValues; Tambient * ones(length(bottomNodes),1)];

    end

    if rightcham == true

        x1r = X_body(5); y1r = Y_body(5);
        x2r = X_body(6); y2r = Y_body(6);

        A = y2r - y1r;
        B = x1r - x2r;
        C = x2r*y1r - x1r*y2r;

        dist = abs(A*X + B*Y + C) / sqrt(A^2 + B^2);

        s = ((X - x1r).*(x2r - x1r) + ...
             (Y - y1r).*(y2r - y1r)) ./ ((x2r - x1r)^2 + (y2r - y1r)^2);

        rightNodes = find(dist < tol & s >= 0 & s <= 1);

        fixedNodes  = [fixedNodes; rightNodes(:)];
        fixedValues = [fixedValues; Tambient * ones(length(rightNodes),1)];

    end

    if topredge == true

        yTop = max(Y_body);
        xMid = (min(X_body) + max(X_body)) / 2;

        topNodes = find(abs(Y - yTop) < tol & X >= xMid);

        fixedNodes  = [fixedNodes; topNodes(:)];
        fixedValues = [fixedValues; Tambient * ones(length(topNodes),1)];

    end


%The normal custom tool

elseif normalMesh == true

    if botedge == true

        bottomNodes = find(abs(Y - min(Yb)) < tol);

        fixedNodes  = [fixedNodes; bottomNodes(:)];
        fixedValues = [fixedValues; Tambient * ones(length(bottomNodes),1)];

    end

    if rightcham == true

        x1r = Xb(3); y1r = Yb(3);
        x2r = Xb(2); y2r = Yb(2);

        A = y2r - y1r;
        B = x1r - x2r;
        C = x2r*y1r - x1r*y2r;

        dist = abs(A*X + B*Y + C) / sqrt(A^2 + B^2);

        s = ((X - x1r).*(x2r - x1r) + ...
             (Y - y1r).*(y2r - y1r)) ./ ((x2r - x1r)^2 + (y2r - y1r)^2);

        rightNodes = find(dist < tol & s >= 0 & s <= 1);

        fixedNodes  = [fixedNodes; rightNodes(:)];
        fixedValues = [fixedValues; Tambient * ones(length(rightNodes),1)];

    end

    if topredge == true

        yTop = max(Yb);
        xMid = (min(Xb) + max(Xb)) / 2;

        topNodes = find(abs(Y - yTop) < tol & X >= xMid);

        fixedNodes  = [fixedNodes; topNodes(:)];
        fixedValues = [fixedValues; Tambient * ones(length(topNodes),1)];

    end

else

    error('No mesh type selected.');

end

%Remove duplicates and clean up

allNodes      = (1:n_nodes)';
ConvNodes     = setdiff(allNodes, fixedNodes);

[fixedNodes, ia] = unique(fixedNodes, 'last');
fixedValues = fixedValues(ia);

end