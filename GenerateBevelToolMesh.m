function [X,Y,ncon,material,n_nodes,n_elements] = GenerateBevelToolMesh( ...
    X_body, Y_body, ...
    X_tri, Y_tri, ...
    nx)

% =========================================================
% GENERATEBEVELTOOLMESH
%
% Geometry (CCW ordering):
%
%   Body polygon (9 pts, last = first):
%     (1) xB1, yB1  = top of upper-left bevel    (0.0002533,  0)
%     (2) xB2, yB2  = bottom of upper-left bevel  (0.0000123, -0.0002121)
%     (3) xB3, yB3  = lower-left bevel bottom     (0.00078,   -0.004)
%     (4) xB4, yB4  = bottom-left                 (0.006,     -0.004)
%     (5) xB5, yB5  = bottom-right                (0.0112,    -0.004)
%     (6) xB6, yB6  = top-right                   (0.012,      0)
%     (7) xB7, yB7  = top inner-right             (0.006,      0)
%     (8) xB8, yB8  = top inner-left              (0.0005533,  0)
%     (9) = (1)  closes
%
%   BUE triangle (4 pts, CCW, last = first):
%     (1) xT1, yT1  = stagnation tip              (0,          0)
%     (2) xT2, yT2  = bottom of bevel             (0.0000123, -0.0002121)
%     (3) xT3, yT3  = top of bevel                (0.0002533,  0)
%     (4) = (1)
%
%   Shared interface:  body edge (1)->(2)  ==  BUE edge (3)->(2)
%
% Mesh strategy:
%   Three independent structured blocks, all sharing nodes on
%   their mutual interfaces via addNode deduplication:
%
%   Block A  (upper body):  quadrilateral region
%              top    y = 0,    x: xB1 -> xB6
%              bottom y = yB2,  x: xB2 -> xRight(yB2)
%              left   bevel segment pt1->pt2  (shared with BUE right edge)
%              right  slant pt6->pt5 (partial)
%
%   Block B  (lower body):  quadrilateral region
%              top    y = yB2,  x: xB2 -> xRight(yB2)
%              bottom y = yB3,  x: xB3 -> xRight(yB3)
%              left   bevel segment pt2->pt3
%              right  slant pt6->pt5 (partial)
%
%   Block C  (BUE triangle):  triangular region
%              base   from T3(xB1,0) -> T2(xB2,yB2)  (= bevel, shared with A)
%              tip    T1(0,0)
%
% material = 1 -> main carbide tool body
% material = 2 -> built-up edge / stagnation zone
% =========================================================

tol = 1e-12;

% ---------------------------------------------------------
% Mesh density
% ---------------------------------------------------------
nx_body  = max(30, 3*nx);   % columns across body width
ny_upper = 8;               % rows in upper body block (y: 0 -> yB2)
ny_lower = 20;              % rows in lower body block (y: yB2 -> yBottom)
nr_bue   = 8;               % rows in BUE triangle (must match ny_upper
                            % so bevel interface nodes are coincident)

% Force nr_bue == ny_upper so the left column of Block A and the
% right edge of the BUE have identical y-levels.
nr_bue   = ny_upper;

% =========================================================
% 1. UNPACK GEOMETRY
% =========================================================
xB1 = X_body(1);  yB1 = Y_body(1);
xB2 = X_body(2);  yB2 = Y_body(2);
xB3 = X_body(3);  yB3 = Y_body(3);
xB4 = X_body(4);  yB4 = Y_body(4);
xB5 = X_body(5);  yB5 = Y_body(5);
xB6 = X_body(6);  yB6 = Y_body(6);

xT1 = X_tri(1);   yT1 = Y_tri(1);   % BUE tip    (0, 0)
xT2 = X_tri(2);   yT2 = Y_tri(2);   % bevel bot  == (xB2, yB2)
xT3 = X_tri(3);   yT3 = Y_tri(3);   % bevel top  == (xB1, yB1)

% =========================================================
% 2. NODE STORE
% =========================================================
X = [];
Y = [];

    function id = addNode(xp, yp)
        existing = find(abs(X - xp) < tol & abs(Y - yp) < tol, 1);
        if isempty(existing)
            X(end+1,1) = xp;
            Y(end+1,1) = yp;
            id = numel(X);
        else
            id = existing;
        end
    end

% =========================================================
% 3. RIGHT BOUNDARY (shared by both body blocks)
% =========================================================
% Single slant: pt6(xB6,yB6=0) -> pt5(xB5,yB5=-0.004)

    function xR = rightBoundaryX(yrow)
        if abs(yB6 - yB5) < tol
            xR = xB6;
        else
            s  = (yrow - yB6) / (yB5 - yB6);
            xR = xB6 + s*(xB5 - xB6);
        end
    end

% =========================================================
% 4. BLOCK A  --  UPPER BODY  (y: 0 -> yB2)
% =========================================================
% Left boundary:  bevel pt1->pt2  (this IS the BUE interface)
% Right boundary: slant pt6->pt5  (partial, from y=0 to y=yB2)
%
% y-levels are shared with the BUE (nr_bue == ny_upper),
% guaranteeing coincident nodes on the bevel edge.

yA = linspace(yB1, yB2, ny_upper);   % top -> bottom, ny_upper levels

nodeA = zeros(ny_upper, nx_body);

for j = 1:ny_upper
    yrow = yA(j);

    % Left x: interpolate along bevel pt1->pt2
    s    = (yrow - yB1) / (yB2 - yB1 + eps);  % eps avoids /0 if equal
    xL   = xB1 + s*(xB2 - xB1);

    xR   = rightBoundaryX(yrow);

    for i = 1:nx_body
        xi = (i-1)/(nx_body-1);
        nodeA(j,i) = addNode((1-xi)*xL + xi*xR, yrow);
    end
end

ncon     = [];
material = [];

for j = 1:(ny_upper-1)
    for i = 1:(nx_body-1)
        n1 = nodeA(j,   i  );
        n2 = nodeA(j,   i+1);
        n3 = nodeA(j+1, i  );
        n4 = nodeA(j+1, i+1);

        ncon     = [ncon;     n1 n2 n3; n2 n4 n3];
        material = [material; 1;        1        ];
    end
end

% =========================================================
% 5. BLOCK B  --  LOWER BODY  (y: yB2 -> yBottom)
% =========================================================
% Left boundary:  bevel pt2->pt3
% Right boundary: slant pt6->pt5  (partial, from y=yB2 to y=yBottom)
% Top row of Block B shares nodes with bottom row of Block A (same y=yB2).

yBottom = yB4;   % -0.004
yB_lev  = linspace(yB2, yBottom, ny_lower);   % ny_lower levels

nodeB = zeros(ny_lower, nx_body);

for j = 1:ny_lower
    yrow = yB_lev(j);

    % Left x: interpolate along bevel pt2->pt3
    s    = (yrow - yB2) / (yB3 - yB2 + eps);
    xL   = xB2 + s*(xB3 - xB2);

    xR   = rightBoundaryX(yrow);

    for i = 1:nx_body
        xi = (i-1)/(nx_body-1);
        nodeB(j,i) = addNode((1-xi)*xL + xi*xR, yrow);
    end
end

for j = 1:(ny_lower-1)
    for i = 1:(nx_body-1)
        n1 = nodeB(j,   i  );
        n2 = nodeB(j,   i+1);
        n3 = nodeB(j+1, i  );
        n4 = nodeB(j+1, i+1);

        ncon     = [ncon;     n1 n2 n3; n2 n4 n3];
        material = [material; 1;        1        ];
    end
end

% =========================================================
% 6. BLOCK C  --  BUE TRIANGLE  (tip T1, base T3->T2)
% =========================================================
% Row 1 (eta=0): base edge T3(xB1,0) -> T2(xB2,yB2)   nr_bue nodes
% Row end (eta=1): tip T1(0,0)                          1 node
%
% Right edge of each BUE row interpolates T2->T1.
% Left  edge of each BUE row interpolates T3->T1.
%
% y-levels match Block A exactly (nr_bue == ny_upper, same linspace),
% so right-edge BUE nodes ARE the left-column Block A nodes.

bueNodeID = cell(nr_bue, 1);

for row = 1:nr_bue
    eta = (row-1)/(nr_bue-1);

    % Left  edge: T3 -> T1
    xRL = (1-eta)*xT3 + eta*xT1;
    yRL = (1-eta)*yT3 + eta*yT1;

    % Right edge: T2 -> T1  (bevel interface, shared with Block A col 1)
    xRR = (1-eta)*xT2 + eta*xT1;
    yRR = (1-eta)*yT2 + eta*yT1;

    nInRow = nr_bue - row + 1;
    bueNodeID{row} = zeros(nInRow, 1);

    for i = 1:nInRow
        xi = 0;
        if nInRow > 1
            xi = (i-1)/(nInRow-1);
        end
        xp = (1-xi)*xRL + xi*xRR;
        yp = (1-xi)*yRL + xi*yRR;
        bueNodeID{row}(i) = addNode(xp, yp);
    end
end

for row = 1:(nr_bue-1)
    upper = bueNodeID{row};
    lower = bueNodeID{row+1};
    nu    = length(upper);
    nl    = length(lower);

    if nl == 1
        for i = 1:(nu-1)
            ncon     = [ncon;     upper(i) upper(i+1) lower(1)];
            material = [material; 2];
        end
    else
        for i = 1:nl
            ncon     = [ncon;     upper(i) upper(i+1) lower(i)];
            material = [material; 2];
            if i < nl
                ncon     = [ncon;     upper(i+1) lower(i+1) lower(i)];
                material = [material; 2];
            end
        end
    end
end

% =========================================================
% 7. MERGE NEAR-DUPLICATE NODES
% =========================================================
coordTol     = 1e-9;
roundedNodes = round([X Y] / coordTol) * coordTol;
[~, uniqueIdx, mapOldToNew] = unique(roundedNodes, 'rows', 'stable');
X    = X(uniqueIdx);
Y    = Y(uniqueIdx);
ncon = mapOldToNew(ncon);

% =========================================================
% 8. REMOVE DEGENERATE TRIANGLES
% =========================================================
Ar = zeros(size(ncon,1), 1);
for e = 1:size(ncon,1)
    x1=X(ncon(e,1)); y1=Y(ncon(e,1));
    x2=X(ncon(e,2)); y2=Y(ncon(e,2));
    x3=X(ncon(e,3)); y3=Y(ncon(e,3));
    Ar(e) = 0.5*abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
end
good     = Ar > 1e-16;
ncon     = ncon(good,:);
material = material(good);

% =========================================================
% 9. REMOVE DUPLICATE TRIANGLES
% =========================================================
[~, ui] = unique(sort(ncon,2), 'rows', 'stable');
ncon     = ncon(ui,:);
material = material(ui);

% =========================================================
% 10. REMOVE UNUSED NODES AND RENUMBER
% =========================================================
usedNodes          = unique(ncon(:));
oldToNew           = zeros(length(X),1);
oldToNew(usedNodes)= 1:length(usedNodes);
X    = X(usedNodes);
Y    = Y(usedNodes);
ncon = oldToNew(ncon);

% =========================================================
% 11. COUNTS AND CONNECTIVITY CHECK
% =========================================================
n_nodes    = length(X);
n_elements = size(ncon,1);

G = graph();
G = addnode(G, n_nodes);
for e = 1:n_elements
    G = addedge(G, ncon(e,1), ncon(e,2));
    G = addedge(G, ncon(e,2), ncon(e,3));
    G = addedge(G, ncon(e,3), ncon(e,1));
end
bins     = conncomp(G);
numParts = max(bins);

fprintf('Connected mesh regions: %d\n', numParts);
fprintf('Total nodes:    %d\n', n_nodes);
fprintf('Total elements: %d\n', n_elements);
fprintf('Tool elements:  %d\n', sum(material==1));
fprintf('BUE elements:   %d\n', sum(material==2));

end


%{
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
%}

