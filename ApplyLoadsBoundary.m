function [F, NDU, dzero, clampedNodes, contactNodes] = ApplyLoadsBoundary( ...
    Ft, Fc, Tc, ...
    botedge, rightcham, topredge, ...
    n_nodes, X, Y, ...
    Xb, Yb, ...
    X_body, Y_body, ...
    X_tri, Y_tri, ...
    bevelled, ANSYSmesh, normalMesh)

% =========================================================
% APPLYLOADSBOUNDARY
%
% Builds:
%   F             global force vector
%   dzero         constrained DOF list
%   NDU           number of constrained DOFs
%   clampedNodes  nodes used for visualising clamps
%   contactNodes  nodes used for visualising load contact
% =========================================================

% ---------------------------------------------------------
% Initialise outputs
% ---------------------------------------------------------
F = zeros(2*n_nodes, 1);
dzero = [];
clampedNodes = [];
contactNodes = [];

tolY = 1e-10;
tol  = 1e-6;

% =========================================================
% 1. CONTACT LOAD NODES
% =========================================================
% Contact line is assumed on top surface y = 0, from x = 0
% to x = contact length Tc.
%{
if bevelled == true
x1 = 0.0002533;  y1 = 0.0;
x2 = 0.0000123;  y2 = -0.0002121;

% Edge vector and length
dx = x2 - x1;
dy = y2 - y1;
L  = sqrt(dx^2 + dy^2);

% For each node, compute:
%   t   = projection parameter along the edge (0=pt1, 1=pt2)
%   d   = perpendicular distance from the line
t = ((X - x1)*dx + (Y - y1)*dy) / L^2;
d = abs((X - x1)*dy - (Y - y1)*dx) / L;

tolD = 1e-9;   % perpendicular distance tolerance (tight — these are mesh nodes)
contactNodes = find(d < tolD & t >= -tolD & t <= 1 + tolD);

else
    %}
contactNodes = find( ...
    X >= 0.0000 & ...
    X <= Tc & ...
    abs(Y - 0.0000) < tolY);

if isempty(contactNodes)
    error('No contact nodes found for contact line load.');
end

% Sort left to right
[~, idx] = sort(X(contactNodes));
contactNodes = contactNodes(idx);

disp('Contact nodes used for line load = ');
disp(contactNodes);

disp('Number of contact nodes = ');
disp(length(contactNodes));

if length(contactNodes) < 2
    error('At least two contact nodes are required for a distributed line load.');
end

% Contact length based on actual selected nodes
xContact = X(contactNodes);
Lc = xContact(end) - xContact(1);

if Lc <= 0
    error('Contact length is zero or negative.');
end

% Force per unit length
qx = Fc / Lc;
qy = Ft / Lc;

% Tributary length for each contact node
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

Fx_nodes = qx * tributary;
Fy_nodes = qy * tributary;

disp('Total Fx applied = ');
disp(sum(Fx_nodes));

disp('Total Fy applied = ');
disp(sum(Fy_nodes));

% Assemble nodal load vector
for i = 1:length(contactNodes)

    n = contactNodes(i);

    F(2*n - 1) = F(2*n - 1) + Fx_nodes(i);
    F(2*n)     = F(2*n)     + Fy_nodes(i);

end


% =========================================================
% 2. BOUNDARY CONDITIONS
% =========================================================

% ---------------------------------------------------------
% ANSYS mesh boundary conditions
% ---------------------------------------------------------
if ANSYSmesh == true

    if botedge == true

        bottomNodes = find(abs(Y - min(Y)) < tol);

        dzero = [dzero;
                 reshape([2*bottomNodes-1; 2*bottomNodes], [], 1)];

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
        den = (x2 - x1)^2 + (y2 - y1)^2;

        s = ((X - x1).*(x2 - x1) + ...
             (Y - y1).*(y2 - y1)) ./ den;

        rightNodes = find(dist < tol & s >= 0 & s <= 1);

        dzero = [dzero;
                 reshape([2*rightNodes-1; 2*rightNodes], [], 1)];

    end

    if topredge == true

        yTop = max(Y);
        xMid = (min(X) + max(X)) / 2;

        topRightNodes = find(abs(Y - yTop) < tol & X >= xMid);

        dzero = [dzero;
                 reshape([2*topRightNodes-1; 2*topRightNodes], [], 1)];

    end


% ---------------------------------------------------------
% Bevelled tool boundary conditions
% ---------------------------------------------------------
elseif bevelled == true

    if botedge == true

        yBottom = min(Y_body);
        bottomNodes = find(abs(Y - yBottom) < tol);

        dzero = [dzero;
                 reshape([2*bottomNodes-1; 2*bottomNodes], [], 1)];

    end

    if rightcham == true

        % Right boundary from top-right to bottom-right
        x1r = X_body(5); y1r = Y_body(5);
        x2r = X_body(6); y2r = Y_body(6);

        A = y2r - y1r;
        B = x1r - x2r;
        C = x2r*y1r - x1r*y2r;

        dist = abs(A*X + B*Y + C) / sqrt(A^2 + B^2);
        den = (x2r - x1r)^2 + (y2r - y1r)^2;

        s = ((X - x1r).*(x2r - x1r) + ...
             (Y - y1r).*(y2r - y1r)) ./ den;

        rightNodes = find(dist < tol & s >= 0 & s <= 1);

        dzero = [dzero;
                 reshape([2*rightNodes-1; 2*rightNodes], [], 1)];

    end

    if topredge == true

        yTop = max(Y_body);
        xMid = (min(X_body) + max(X_body)) / 2;

        topRightNodes = find(abs(Y - yTop) < tol & X >= xMid);

        dzero = [dzero;
                 reshape([2*topRightNodes-1; 2*topRightNodes], [], 1)];

    end


% ---------------------------------------------------------
% Normal custom tool boundary conditions
% ---------------------------------------------------------
elseif normalMesh == true

    if botedge == true

        yBottom = min(Yb);
        bottomNodes = find(abs(Y - yBottom) < tol);

        dzero = [dzero;
                 reshape([2*bottomNodes-1; 2*bottomNodes], [], 1)];

    end

    if rightcham == true

        x1r = Xb(3); y1r = Yb(3);
        x2r = Xb(2); y2r = Yb(2);

        A = y2r - y1r;
        B = x1r - x2r;
        C = x2r*y1r - x1r*y2r;

        dist = abs(A*X + B*Y + C) / sqrt(A^2 + B^2);
        den = (x2r - x1r)^2 + (y2r - y1r)^2;

        s = ((X - x1r).*(x2r - x1r) + ...
             (Y - y1r).*(y2r - y1r)) ./ den;

        rightNodes = find(dist < tol & s >= 0 & s <= 1);

        dzero = [dzero;
                 reshape([2*rightNodes-1; 2*rightNodes], [], 1)];

    end

    if topredge == true

        yTop = max(Yb);
        xMid = (min(Xb) + max(Xb)) / 2;

        topRightNodes = find(abs(Y - yTop) < tol & X >= xMid);

        dzero = [dzero;
                 reshape([2*topRightNodes-1; 2*topRightNodes], [], 1)];

    end

else

    error('No mesh type selected: ANSYSmesh, bevelled, and normalMesh are all false.');

end


% =========================================================
% 3. FINALISE CONSTRAINT LIST
% =========================================================
dzero = unique(dzero);
NDU = length(dzero);
clampedNodes = unique(ceil(dzero/2));

fprintf('Number of constrained DOFs NDU = %d\n', NDU);
fprintf('Number of clamped nodes = %d\n', length(clampedNodes));

if isempty(dzero)
    error('No constrained DOFs found. Boundary conditions were not applied.');
end

end