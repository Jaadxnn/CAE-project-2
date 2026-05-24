function [X,Y,ncon,material,n_nodes,n_elements] = GenerateBevelToolMesh( ...
    X_body, Y_body, ...
    X_tri, Y_tri, ...
    nx)

% =========================================================
% GENERATEBEVELTOOLMESH
%
% Structured/semi-structured triangular mesh for bevelled
% cutting tool geometry.
%
% material = 1 -> main carbide tool body
% material = 2 -> built-up edge / stagnation zone
% =========================================================

tol = 1e-12;

% ---------------------------------------------------------
% Mesh density settings
% ---------------------------------------------------------
% The app currently passes nx = 10, which is too coarse.
% We use it as a minimum guide but enforce a better mesh.
nx_body = max(30, 3*nx);
ny_body = 16;

nr_bue = 8;     % number of BUE rows, including top and bottom

% =========================================================
% 1. IDENTIFY BODY GEOMETRY POINTS
% =========================================================
% Current bevelled body coordinates:
%
% X_body(1),Y_body(1) = top-left bevel point
% X_body(2),Y_body(2) = top-right
% X_body(3),Y_body(3) = bottom-right
% X_body(4),Y_body(4) = bottom-left
% X_body(5),Y_body(5) = lower-left bevel point
% X_body(6),Y_body(6) = close point
%
% Current BUE triangle:
%
% X_tri(1),Y_tri(1) = sharp stagnation tip
% X_tri(2),Y_tri(2) = top bevel point
% X_tri(3),Y_tri(3) = lower bevel point
% X_tri(4),Y_tri(4) = close point

xTopLeft  = X_body(1);
yTopLeft  = Y_body(1);

xTopRight = X_body(2);
yTopRight = Y_body(2);

xBotRight = X_body(3);
yBotRight = Y_body(3);

xBotLeft  = X_body(4);
yBotLeft  = Y_body(4);

xChamLow  = X_body(5);
yChamLow  = Y_body(5);

% BUE triangle points
xA = X_tri(1); yA = Y_tri(1);   % tip
xB = X_tri(2); yB = Y_tri(2);   % top bevel point
xC = X_tri(3); yC = Y_tri(3);   % lower bevel point

% =========================================================
% 2. CREATE BODY Y-LEVELS
% =========================================================
% Include BUE y-levels so the BUE and body share nodes along
% the bevel interface.
% =========================================================

yMain = linspace(yTopLeft, yBotLeft, ny_body);
yBUE  = linspace(yB, yC, nr_bue);

yLevels = unique([yMain yBUE], 'stable');
yLevels = sort(yLevels, 'descend');

% =========================================================
% 3. HELPER: FIND OR ADD NODE
% =========================================================
X = [];
Y = [];

    function id = addNode(xp, yp)
        existing = find(abs(X - xp) < tol & abs(Y - yp) < tol, 1);

        if isempty(existing)
            X = [X; xp];
            Y = [Y; yp];
            id = length(X);
        else
            id = existing;
        end
    end

% =========================================================
% 4. MAIN BODY MESH
% =========================================================
% The left boundary is piecewise:
%   top-left bevel point -> lower-left bevel point
%   lower-left bevel point -> bottom-left
%
% The right boundary is:
%   top-right -> bottom-right
% =========================================================

bodyNodeID = zeros(length(yLevels), nx_body);

for j = 1:length(yLevels)

    yrow = yLevels(j);

    % -------------------------------
    % Left boundary x-value
    % -------------------------------
    if yrow >= yChamLow
        % On bevel edge from top-left to lower-left bevel point
        if abs(yTopLeft - yChamLow) < tol
            sL = 0;
        else
            sL = (yrow - yTopLeft) / (yChamLow - yTopLeft);
        end

        xLeft = xTopLeft + sL*(xChamLow - xTopLeft);
        yLeft = yrow;

    else
        % On left lower edge from lower-left bevel point to bottom-left
        if abs(yChamLow - yBotLeft) < tol
            sL = 0;
        else
            sL = (yrow - yChamLow) / (yBotLeft - yChamLow);
        end

        xLeft = xChamLow + sL*(xBotLeft - xChamLow);
        yLeft = yrow;
    end

    % -------------------------------
    % Right boundary x-value
    % -------------------------------
    if abs(yTopRight - yBotRight) < tol
        sR = 0;
    else
        sR = (yrow - yTopRight) / (yBotRight - yTopRight);
    end

    xRight = xTopRight + sR*(xBotRight - xTopRight);
    yRight = yrow;

    % -------------------------------
    % Fill row from left boundary to right boundary
    % -------------------------------
    for i = 1:nx_body

        xi = (i-1)/(nx_body-1);

        xp = (1-xi)*xLeft + xi*xRight;
        yp = (1-xi)*yLeft + xi*yRight;

        bodyNodeID(j,i) = addNode(xp, yp);

    end
end

% Body connectivity
ncon = [];
material = [];

for j = 1:(size(bodyNodeID,1)-1)

    for i = 1:(nx_body-1)

        n1 = bodyNodeID(j,i);
        n2 = bodyNodeID(j,i+1);
        n3 = bodyNodeID(j+1,i);
        n4 = bodyNodeID(j+1,i+1);

        % Split each quad into two triangles
        ncon = [ncon; n1 n2 n3];
        material = [material; 1];

        ncon = [ncon; n2 n4 n3];
        material = [material; 1];

    end
end

% =========================================================
% 5. BUE TRIANGLE MESH
% =========================================================
% Rows go from top edge A-B down to point C.
% The right edge B-C is shared with the tool body bevel edge.
% =========================================================

bueNodeID = cell(nr_bue,1);

for row = 1:nr_bue

    eta = (row-1)/(nr_bue-1);

    % Left edge A -> C
    xLeft = (1-eta)*xA + eta*xC;
    yLeft = (1-eta)*yA + eta*yC;

    % Right edge B -> C
    xRight = (1-eta)*xB + eta*xC;
    yRight = (1-eta)*yB + eta*yC;

    % Number of nodes decreases towards the triangle tip
    nInRow = nr_bue - row + 1;

    bueNodeID{row} = zeros(nInRow,1);

    for i = 1:nInRow

        if nInRow == 1
            xi = 0;
        else
            xi = (i-1)/(nInRow-1);
        end

        xp = (1-xi)*xLeft + xi*xRight;
        yp = (1-xi)*yLeft + xi*yRight;

        bueNodeID{row}(i) = addNode(xp, yp);

    end
end

% BUE connectivity
for row = 1:(nr_bue-1)

    upper = bueNodeID{row};
    lower = bueNodeID{row+1};

    nu = length(upper);
    nl = length(lower);

    % Since lower row has one fewer node than upper row
    for i = 1:nl

        % Upper triangle
        n1 = upper(i);
        n2 = upper(i+1);
        n3 = lower(i);

        ncon = [ncon; n1 n2 n3];
        material = [material; 2];

        % Lower connecting triangle, where possible
        if i < nl
            n4 = lower(i+1);

            ncon = [ncon; n2 n4 n3];
            material = [material; 2];
        end
    end
end

% =========================================================
% 6. MERGE NEAR-DUPLICATE NODES
% =========================================================
% The BUE and tool body may create nodes that are visually
% coincident but numerically slightly different. If they are
% not merged, the BUE can become disconnected from the tool
% body, causing a singular stiffness matrix.

coordTol = 1e-9;

roundedNodes = round([X Y] / coordTol) * coordTol;

[~, uniqueIdx, mapOldToNew] = unique(roundedNodes, 'rows', 'stable');

X = X(uniqueIdx);
Y = Y(uniqueIdx);

ncon = mapOldToNew(ncon);


% =========================================================
% 7. REMOVE DEGENERATE TRIANGLES
% =========================================================
A = zeros(size(ncon,1),1);

for e = 1:size(ncon,1)

    x1 = X(ncon(e,1)); y1 = Y(ncon(e,1));
    x2 = X(ncon(e,2)); y2 = Y(ncon(e,2));
    x3 = X(ncon(e,3)); y3 = Y(ncon(e,3));

    A(e) = 0.5 * abs( ...
        x1*(y2-y3) + ...
        x2*(y3-y1) + ...
        x3*(y1-y2));

end

good = A > 1e-16;

ncon = ncon(good,:);
material = material(good);


% =========================================================
% 8. REMOVE DUPLICATE TRIANGLES
% =========================================================
% Duplicate triangles can occur after node merging.

sortedCon = sort(ncon,2);

[~, uniqueElemIdx] = unique(sortedCon, 'rows', 'stable');

ncon = ncon(uniqueElemIdx,:);
material = material(uniqueElemIdx);


% =========================================================
% 9. REMOVE UNUSED NODES AND RENUMBER CONNECTIVITY
% =========================================================
usedNodes = unique(ncon(:));

oldToNew = zeros(length(X),1);
oldToNew(usedNodes) = 1:length(usedNodes);

X = X(usedNodes);
Y = Y(usedNodes);

ncon = oldToNew(ncon);


% =========================================================
% 10. FINAL OUTPUT COUNTS
% =========================================================
n_nodes = length(X);
n_elements = size(ncon,1);

% =========================================================
% 11. BASIC CONNECTIVITY CHECK
% =========================================================
G = graph();

G = addnode(G, n_nodes);

for e = 1:size(ncon,1)
    G = addedge(G, ncon(e,1), ncon(e,2));
    G = addedge(G, ncon(e,2), ncon(e,3));
    G = addedge(G, ncon(e,3), ncon(e,1));
end

bins = conncomp(G);
numParts = max(bins);

fprintf('Connected mesh regions: %d\n', numParts);

fprintf('Total nodes:    %d\n', n_nodes);
fprintf('Total elements: %d\n', n_elements);
fprintf('Tool elements:  %d\n', sum(material == 1));
fprintf('BUE elements:   %d\n', sum(material == 2));

end