function boundaryNodes = findBoundaryNodes(ncon)

% =========================================================
% FINDBOUNDARYNODES
%
% Finds external boundary nodes from triangular connectivity.
% Boundary edges are edges that appear only once.
% =========================================================

edges = [ ...
    ncon(:,[1 2]);
    ncon(:,[2 3]);
    ncon(:,[3 1])];

% Sort node numbers in each edge so [2 5] and [5 2] match
edges = sort(edges, 2);

% Count repeated edges
[uniqueEdges, ~, ic] = unique(edges, 'rows');
edgeCount = accumarray(ic, 1);

% Boundary edges appear only once
boundaryEdges = uniqueEdges(edgeCount == 1, :);

% Boundary nodes are all nodes in boundary edges
boundaryNodes = unique(boundaryEdges(:));

end