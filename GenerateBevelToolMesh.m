function [X,Y,ncon,material,n_nodes,n_elements] = GenerateBevelToolMesh( ...
    X_body, Y_body, ...
    X_tri, Y_tri, ...
    nx)


tol = 1e-12;

nx_body = max(30, 3*nx);%columns across body width
ny_upper = 8;%rows in upper body block
ny_lower = 20;%rows in lower body block (y: yB2 -> yBottom)
nr_bue = 8;%rows in BUE triangle (must match ny_upper so bevel interface nodes coincide)

%Force nr_bue to match ny_upper so the left column of Block A and the
%right edge of the BUE have identical y levels.
nr_bue   = ny_upper;

%Geometrey of the bevelled tool
xB1 = X_body(1);  yB1 = Y_body(1);
xB2 = X_body(2);  yB2 = Y_body(2);
xB3 = X_body(3);  yB3 = Y_body(3);
xB4 = X_body(4);  yB4 = Y_body(4);
xB5 = X_body(5);  yB5 = Y_body(5);
xB6 = X_body(6);  yB6 = Y_body(6);

xT1 = X_tri(1);   yT1 = Y_tri(1);%BUE tip (0, 0)
xT2 = X_tri(2);   yT2 = Y_tri(2);%bevel bot==(xB2, yB2)
xT3 = X_tri(3);   yT3 = Y_tri(3);%bevel top ==(xB1, yB1)

%Store the nodes
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


%BLOCK A upper tool region
% Structured mesh between the bevel edge and the right boundary.
% This region shares its left boundary with the BUE interface.

yA = linspace(yB1, yB2, ny_upper); %Top to bottom, ny_upper levels

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

%Connectivty, so split the quads into there triagnles
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

%BLOCK B the lower tool region
%Same structured approach as Block A but below bevel.
%Shares interface nodes with upper region for continuity.

yBottom = yB4;%-0.004
yB_lev  = linspace(yB2, yBottom, ny_lower); 

nodeB = zeros(ny_lower, nx_body);

for j = 1:ny_lower
    yrow = yB_lev(j);

    %Left x interpolate along bevel pt2 to pt3
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

%BLOCK C: This is the BUILT UP EDGE (BUE) triangle
%Row 1 (eta=0): base edge T3(xB1,0) -> T2(xB2,yB2) nr_bue nodes
%Row end (eta=1): tip T1(0,0) 1 node
%Right edge of each BUE row interpolates T2->T1.
%Left  edge of each BUE row interpolates T3->T1.
%y-levels match Block A exactly (nr_bue == ny_upper, same linspace),
%so right-edge BUE nodes ARE the left-column Block A nodes.

bueNodeID = cell(nr_bue, 1);

for row = 1:nr_bue
    eta = (row-1)/(nr_bue-1);

    % Left and right edges of the BUE
    xRL = (1-eta)*xT3 + eta*xT1;
    yRL = (1-eta)*yT3 + eta*yT1;


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

%Merge nodes that are pretty much duplicates. So any nodes that are 1e-9
%apart are considered at the same point. E.g node 47 and 317 are considered
%duplicates they are merged to node 47 and the elemnt that once pointed
%to317 now point to 47 which acts as a stitching mechanism between the BUE
%and the main body
coordTol     = 1e-9;
roundedNodes = round([X Y] / coordTol) * coordTol;
[~, uniqueIdx, mapOldToNew] = unique(roundedNodes, 'rows', 'stable');
X    = X(uniqueIdx);
Y    = Y(uniqueIdx);
ncon = mapOldToNew(ncon);

%Remove degenerate traingles. triangles that have no area
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

%Remove  duplicated triangles
[~, ui] = unique(sort(ncon,2), 'rows', 'stable');
ncon     = ncon(ui,:);
material = material(ui);


%Remove duplicate or unsued nodes then shuffle
usedNodes          = unique(ncon(:));
oldToNew           = zeros(length(X),1);
oldToNew(usedNodes)= 1:length(usedNodes);
X    = X(usedNodes);
Y    = Y(usedNodes);
ncon = oldToNew(ncon);


n_nodes    = length(X);
n_elements = size(ncon,1);

fprintf('Total nodes:    %d\n', n_nodes);
fprintf('Total elements: %d\n', n_elements);
fprintf('Tool elements:  %d\n', sum(material==1));
fprintf('BUE elements:   %d\n', sum(material==2));

end



