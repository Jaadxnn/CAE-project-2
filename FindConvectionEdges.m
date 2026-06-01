function convEdges = FindConvectionEdges(ncon, X, Y, segments)

tol = 1e-6;

% ---------------------------------------------------------
% Step 1: Find all boundary edges
% An edge is on the boundary if it belongs to only 1 element
% ---------------------------------------------------------

% List every edge of every element
allEdges = [ncon(:,[1,2]); ncon(:,[2,3]); ncon(:,[3,1])];

% Sort each edge so [i,j] and [j,i] are treated the same
allEdges = sort(allEdges, 2);

% Count how many times each unique edge appears
[uniqueEdges, ~, ic] = unique(allEdges, 'rows');
counts               = accumarray(ic, 1);

% Keep only edges that appear once (boundary edges)
bndEdges = uniqueEdges(counts == 1, :);

% ---------------------------------------------------------
% Step 2: Check which boundary edges lie on our segments
% ---------------------------------------------------------
convEdges = [];

for s = 1:size(segments, 1)

    x1 = segments(s,1);  y1 = segments(s,2);
    x2 = segments(s,3);  y2 = segments(s,4);
    dx = x2 - x1;
    dy = y2 - y1;

    % For each node, check if it lies on this segment
    dist  = abs(dy*X - dx*Y + x2*y1 - y2*x1) / sqrt(dx^2 + dy^2);
    t     = ((X - x1)*dx + (Y - y1)*dy)       / (dx^2 + dy^2);
    onSeg = dist < tol & t >= 0 & t <= 1;

    % Keep boundary edges where both nodes are on the segment
    for e = 1:size(bndEdges, 1)
        ni = bndEdges(e,1);
        nj = bndEdges(e,2);
        if onSeg(ni) && onSeg(nj)
            convEdges = [convEdges; ni nj];
        end
    end

end

end