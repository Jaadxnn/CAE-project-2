function [X,Y,ncon,material, n_nodes, n_elements] = GenerateBevelToolMesh( ...
    X_body, Y_body, ...
    X_tri, Y_tri, ...
    nx)

%% ---------------------------------------------------------
% DOMAIN LIMITS
% ---------------------------------------------------------
xmin = min([X_body; X_tri]);
xmax = max([X_body; X_tri]);
ymin = min([Y_body; Y_tri]);
ymax = max([Y_body; Y_tri]);

%% ---------------------------------------------------------
% HEXAGONAL POINT LATTICE
% ---------------------------------------------------------
dx = (xmax - xmin) / (nx - 1);
dy = dx * sqrt(3) / 2;

y_rows = ymin : dy : ymax;
xg = [];
yg = [];

for row = 1:length(y_rows)
    if mod(row, 2) == 0
        x_vals = (xmin + dx/2) : dx : xmax;
    else
        x_vals = xmin : dx : xmax;
    end
    xg = [xg; x_vals(:)];
    yg = [yg; repmat(y_rows(row), numel(x_vals), 1)];
end

%% ---------------------------------------------------------
% KEEP ONLY INTERIOR POINTS INSIDE EITHER REGION
% ---------------------------------------------------------
inside_body = inpolygon(xg, yg, X_body, Y_body);
inside_tri  = inpolygon(xg, yg, X_tri,  Y_tri);
inside      = inside_body | inside_tri;

xg = xg(inside);
yg = yg(inside);

%% ---------------------------------------------------------
% SEED INTERIOR POINTS INSIDE BUE TRIANGLE
% ---------------------------------------------------------
n_bue_seeds = 8;
bue_pts_x = zeros(n_bue_seeds, 1);
bue_pts_y = zeros(n_bue_seeds, 1);

vx = X_tri(1:end-1);
vy = Y_tri(1:end-1);

k = 0;
for i = 1:4
    for j = 1:4
        a1 = i/5;
        a2 = j/5;
        a3 = 1 - a1 - a2;
        if a3 > 0
            k = k + 1;
            bue_pts_x(k) = a1*vx(1) + a2*vx(2) + a3*vx(3);
            bue_pts_y(k) = a1*vy(1) + a2*vy(2) + a3*vy(3);
        end
    end
end
bue_pts_x = bue_pts_x(1:k);
bue_pts_y = bue_pts_y(1:k);

%% ---------------------------------------------------------
% COMBINE BOUNDARY + INTERIOR NODES
% ---------------------------------------------------------
X = [X_body; X_tri(1:end-1); xg; bue_pts_x];
Y = [Y_body; Y_tri(1:end-1); yg; bue_pts_y];

%% ---------------------------------------------------------
% REMOVE DUPLICATE NODES
% ---------------------------------------------------------
nodes = unique([X Y], 'rows');
X = nodes(:,1);
Y = nodes(:,2);

%% ---------------------------------------------------------
% DELAUNAY TRIANGULATION
% ---------------------------------------------------------
DT   = delaunayTriangulation(X, Y);
ncon = DT.ConnectivityList;

%% ---------------------------------------------------------
% REMOVE TRIANGLES OUTSIDE BOTH REGIONS
% ---------------------------------------------------------
xc = mean(X(ncon), 2);
yc = mean(Y(ncon), 2);

inside_body = inpolygon(xc, yc, X_body, Y_body);
inside_tri  = inpolygon(xc, yc, X_tri,  Y_tri);
inside      = inside_body | inside_tri;

ncon = ncon(inside, :);

%% ---------------------------------------------------------
% ASSIGN MATERIALS BY ELEMENT CENTROID
% ---------------------------------------------------------
xc = mean(X(ncon), 2);
yc = mean(Y(ncon), 2);

is_tri   = inpolygon(xc, yc, X_tri, Y_tri);
material = ones(size(ncon,1), 1);
material(is_tri) = 2;

%% ---------------------------------------------------------
% REMOVE DEGENERATE TRIANGLES
% ---------------------------------------------------------
A = zeros(size(ncon,1), 1);
for e = 1:size(ncon,1)
    x1 = X(ncon(e,1));  y1 = Y(ncon(e,1));
    x2 = X(ncon(e,2));  y2 = Y(ncon(e,2));
    x3 = X(ncon(e,3));  y3 = Y(ncon(e,3));
    A(e) = 0.5 * abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
end

good     = A > 1e-14;
ncon     = ncon(good, :);
material = material(good);

n_nodes    = length(X);
n_elements = size(ncon, 1);

%% ---------------------------------------------------------
% VERIFICATION
% ---------------------------------------------------------
fprintf('Total nodes:    %d\n', n_nodes);
fprintf('Total elements: %d\n', n_elements);
fprintf('Tool elements:  %d\n', sum(material == 1));
fprintf('BUE elements:   %d\n', sum(material == 2));

end