% function [X, Y, ncon, material, n_nodes, n_elements] = GenerateBevelToolMesh( ...
%     Xb, Yb, Xt, Yt, L_ref, nx_ref, nx_coarse, ny, x_refine_start, x_refine_end)
% 
% % GenerateBevelToolMesh — unstructured Delaunay triangular mesh
% % Outputs: X, Y, ncon, material, n_nodes, n_elements (unchanged signature)
% 
% %% 1. SETUP
% Xb = Xb(:);  Yb = Yb(:);
% Xt = Xt(:);  Yt = Yt(:);
% nx_ref    = max(6,  round(nx_ref));
% nx_coarse = max(4,  round(nx_coarse));
% 
% xmin = min([Xb; Xt]);  xmax = max([Xb; Xt]);
% ymin = min([Yb; Yt]);  ymax = max([Yb; Yt]);
% W = xmax - xmin;
% H = ymax - ymin;
% 
% % Fine zone bounds
% if nargin >= 10 && ~isempty(x_refine_start) && ~isempty(x_refine_end) ...
%         && ~any(isnan([x_refine_start, x_refine_end]))
%     xf0 = max(xmin, x_refine_start);
%     xf1 = min(xmax, x_refine_end);
% else
%     xf0 = xmin;
%     xf1 = min(xmin + L_ref, xmax);
% end
% 
% h_fine   = (xf1 - xf0) / nx_ref;
% h_coarse = max((xmax - xmin) / (nx_ref + nx_coarse), h_fine * 1.5);
% 
% fprintf('GenerateBevelToolMesh: fine=[%.4g, %.4g]  h_fine=%.3e  h_coarse=%.3e\n', ...
%         xf0, xf1, h_fine, h_coarse);
% 
% %% 2. LOCAL SPACING  h(x)  — piecewise linear transition
% tw = 3 * h_coarse;  % transition width
% 
% function hv = hx(xi)
%     if xi >= xf0 && xi <= xf1
%         hv = h_fine;
%     elseif xi < xf0
%         hv = h_fine + (h_coarse - h_fine) * min(1, (xf0 - xi) / tw);
%     else
%         hv = h_fine + (h_coarse - h_fine) * min(1, (xi - xf1) / tw);
%     end
% end
% 
% %% 3. BOUNDARY NODES — sample every polygon edge at local h(x)
% X_bnd = [];  Y_bnd = [];
% for poly = 1:2
%     if poly == 1;  Xp = Xb;  Yp = Yb;
%     else;          Xp = Xt;  Yp = Yt;
%     end
%     nv = numel(Xp);
%     for ii = 1:nv
%         jj   = mod(ii, nv) + 1;
%         elen = hypot(Xp(jj)-Xp(ii), Yp(jj)-Yp(ii));
%         npt  = max(3, round(elen / hx(0.5*(Xp(ii)+Xp(jj)))));
%         tt   = linspace(0,1,npt)';
%         X_bnd = [X_bnd;  Xp(ii) + tt*(Xp(jj)-Xp(ii))];  %#ok<AGROW>
%         Y_bnd = [Y_bnd;  Yp(ii) + tt*(Yp(jj)-Yp(ii))];  %#ok<AGROW>
%     end
% end
% 
% %% 4. INTERIOR SEEDS — triangular lattice with jitter in x AND y
% %    Triangular lattice (offset alternate rows by h/2) breaks the
% %    column structure that made the previous mesh look structured.
% rng(42);
% jitter = 0.35;   % fraction of h — high enough to randomise, low enough to stay dense
% 
% X_int = [];  Y_int = [];
% row   = 0;
% y_cur = ymin;
% 
% while y_cur <= ymax + h_coarse
%     % x-walk for this row with graded spacing
%     x_off = mod(row, 2) * 0.5;   % stagger odd rows by half a pitch
%     x_cur = xmin + x_off * hx(xmin);
% 
%     while x_cur <= xmax + h_fine
%         hc = hx(x_cur);
%         % Jitter in BOTH x and y
%         xj = x_cur + (rand - 0.5) * jitter * hc;
%         yj = y_cur + (rand - 0.5) * jitter * hc;
%         X_int(end+1,1) = xj;  %#ok<AGROW>
%         Y_int(end+1,1) = yj;  %#ok<AGROW>
%         x_cur = x_cur + hc;
%     end
% 
%     % Vertical step: sqrt(3)/2 * h for triangular lattice
%     h_mid = hx(xmin + W/2);
%     y_cur = y_cur + h_mid * 0.866;
%     row   = row + 1;
% end
% 
% %% 5. ASSEMBLE, FILTER, DEDUPLICATE
% X_all = [X_bnd; Xb; Xt; X_int];
% Y_all = [Y_bnd; Yb; Yt; Y_int];
% 
% % Strict polygon filter — removes anything outside domain
% in_b  = inpolygon(X_all, Y_all, Xb, Yb);
% in_t  = inpolygon(X_all, Y_all, Xt, Yt);
% X_all = [X_all(in_b | in_t); Xb; Xt];
% Y_all = [Y_all(in_b | in_t); Yb; Yt];
% 
% sc  = max(W, H);
% tol = sc * 1e-9;
% [~, ia] = unique(round([X_all, Y_all] / tol) * tol, 'rows', 'stable');
% X_all = X_all(ia);
% Y_all = Y_all(ia);
% 
% %% 6. DELAUNAY TRIANGULATION
% DT   = delaunayTriangulation(X_all, Y_all);
% ncon = DT.ConnectivityList;
% X    = DT.Points(:,1);
% Y    = DT.Points(:,2);
% 
% %% 7. REMOVE EXTERIOR & DEGENERATE TRIANGLES
% xce  = (X(ncon(:,1)) + X(ncon(:,2)) + X(ncon(:,3))) / 3;
% yce  = (Y(ncon(:,1)) + Y(ncon(:,2)) + Y(ncon(:,3))) / 3;
% 
% in_b = inpolygon(xce, yce, Xb, Yb);
% in_t = inpolygon(xce, yce, Xt, Yt);
% keep = in_b | in_t;
% ncon = ncon(keep,:);
% xce  = xce(keep);
% yce  = yce(keep);
% 
% A    = 0.5 * abs(X(ncon(:,1)).*(Y(ncon(:,2))-Y(ncon(:,3))) + ...
%                  X(ncon(:,2)).*(Y(ncon(:,3))-Y(ncon(:,1))) + ...
%                  X(ncon(:,3)).*(Y(ncon(:,1))-Y(ncon(:,2))));
% good = A > 1e-14 * sc^2;
% ncon = ncon(good,:);
% xce  = xce(good);
% yce  = yce(good);
% 
% %% 8. REINDEX
% used      = unique(ncon(:));
% map       = zeros(numel(X),1);
% map(used) = 1:numel(used);
% X    = X(used);
% Y    = Y(used);
% ncon = map(ncon);
% 
% %% 9. MATERIAL  (1=body, 2=insert)
% material = ones(size(ncon,1),1);
% material(inpolygon(xce, yce, Xt, Yt)) = 2;
% 
% %% 10. REPORT
% n_nodes    = numel(X);
% n_elements = size(ncon,1);
% fprintf('GenerateBevelToolMesh: %d nodes, %d elements\n', n_nodes, n_elements);
% fprintf('  Material 1 (body): %d   Material 2 (insert): %d\n', ...
%         sum(material==1), sum(material==2));
% end



function [X, Y, ncon, material, n_nodes, n_elements] = GenerateBevelToolMesh2( ...
    Xb, Yb, Xt, Yt, L_ref, nx_ref, nx_coarse, ny, x_refine_start, x_refine_end)

% GenerateBevelToolMesh — unstructured Delaunay triangular mesh
% Outputs: X, Y, ncon, material, n_nodes, n_elements (unchanged signature)

%% 1. SETUP
Xb = Xb(:);  Yb = Yb(:);
Xt = Xt(:);  Yt = Yt(:);
nx_ref    = max(6,  round(nx_ref));
nx_coarse = max(4,  round(nx_coarse));

xmin = min([Xb; Xt]);  xmax = max([Xb; Xt]);
ymin = min([Yb; Yt]);  ymax = max([Yb; Yt]);
W = xmax - xmin;
H = ymax - ymin;

% Fine zone bounds
if nargin >= 10 && ~isempty(x_refine_start) && ~isempty(x_refine_end) ...
        && ~any(isnan([x_refine_start, x_refine_end]))
    xf0 = max(xmin, x_refine_start);
    xf1 = min(xmax, x_refine_end);
else
    xf0 = xmin;
    xf1 = min(xmin + L_ref, xmax);
end

h_fine   = (xf1 - xf0) / nx_ref;
h_coarse = max((xmax - xmin) / (nx_ref + nx_coarse), h_fine * 1.5);

fprintf('GenerateBevelToolMesh: fine=[%.4g, %.4g]  h_fine=%.3e  h_coarse=%.3e\n', ...
        xf0, xf1, h_fine, h_coarse);

%% 2. LOCAL SPACING  h(x)  — piecewise linear transition
tw = 5 * h_coarse;  % transition half-width — wider = smoother seam

function hv = hx(xi)
    % Cosine blend: smooth S-curve between h_fine and h_coarse.
    % Inside fine zone → h_fine.  Outside transition → h_coarse.
    % In between → smooth cosine ramp (no kink at the boundary).
    if xi >= xf0 && xi <= xf1
        hv = h_fine;
    elseif xi < xf0
        t  = min(1, (xf0 - xi) / tw);          % 0 at edge, 1 far away
        hv = h_fine + (h_coarse - h_fine) * 0.5 * (1 - cos(pi * t));
    else
        t  = min(1, (xi - xf1) / tw);
        hv = h_fine + (h_coarse - h_fine) * 0.5 * (1 - cos(pi * t));
    end
end

%% 3. BOUNDARY NODES — sample every polygon edge at local h(x)
X_bnd = [];  Y_bnd = [];
for poly = 1:2
    if poly == 1;  Xp = Xb;  Yp = Yb;
    else;          Xp = Xt;  Yp = Yt;
    end
    nv = numel(Xp);
    for ii = 1:nv
        jj   = mod(ii, nv) + 1;
        elen = hypot(Xp(jj)-Xp(ii), Yp(jj)-Yp(ii));
        npt  = max(3, round(elen / hx(0.5*(Xp(ii)+Xp(jj)))));
        tt   = linspace(0,1,npt)';
        X_bnd = [X_bnd;  Xp(ii) + tt*(Xp(jj)-Xp(ii))];  %#ok<AGROW>
        Y_bnd = [Y_bnd;  Yp(ii) + tt*(Yp(jj)-Yp(ii))];  %#ok<AGROW>
    end
end

%% 4. INTERIOR SEEDS — triangular lattice with jitter in x AND y
%    Triangular lattice (offset alternate rows by h/2) breaks the
%    column structure that made the previous mesh look structured.
rng(42);
jitter = 0.35;   % fraction of h — high enough to randomise, low enough to stay dense

X_int = [];  Y_int = [];
row   = 0;
y_cur = ymin;

while y_cur <= ymax + h_coarse
    x_off = mod(row, 2) * 0.5;
    x_cur = xmin + x_off * hx(xmin);

    h_row = 0;  n_row = 0;
    while x_cur <= xmax + h_fine
        hc    = hx(x_cur);
        h_row = h_row + hc;
        n_row = n_row + 1;
        xj = x_cur + (rand - 0.5) * jitter * hc;
        yj = y_cur + (rand - 0.5) * jitter * hc;
        X_int(end+1,1) = xj;  %#ok<AGROW>
        Y_int(end+1,1) = yj;  %#ok<AGROW>
        x_cur = x_cur + hc;
    end

    % y-step: average h across this row keeps vertical density consistent
    h_avg = h_row / max(n_row, 1);
    y_cur = y_cur + h_avg * 0.866;
    row   = row + 1;
end

%% 5. ASSEMBLE, FILTER, DEDUPLICATE
X_all = [X_bnd; Xb; Xt; X_int];
Y_all = [Y_bnd; Yb; Yt; Y_int];

% Strict polygon filter — removes anything outside domain
in_b  = inpolygon(X_all, Y_all, Xb, Yb);
in_t  = inpolygon(X_all, Y_all, Xt, Yt);
X_all = [X_all(in_b | in_t); Xb; Xt];
Y_all = [Y_all(in_b | in_t); Yb; Yt];

sc  = max(W, H);
tol = sc * 1e-9;
[~, ia] = unique(round([X_all, Y_all] / tol) * tol, 'rows', 'stable');
X_all = X_all(ia);
Y_all = Y_all(ia);

%% 6. DELAUNAY TRIANGULATION
DT   = delaunayTriangulation(X_all, Y_all);
ncon = DT.ConnectivityList;
X    = DT.Points(:,1);
Y    = DT.Points(:,2);

%% 7. REMOVE EXTERIOR & DEGENERATE TRIANGLES
xce  = (X(ncon(:,1)) + X(ncon(:,2)) + X(ncon(:,3))) / 3;
yce  = (Y(ncon(:,1)) + Y(ncon(:,2)) + Y(ncon(:,3))) / 3;

in_b = inpolygon(xce, yce, Xb, Yb);
in_t = inpolygon(xce, yce, Xt, Yt);
keep = in_b | in_t;
ncon = ncon(keep,:);
xce  = xce(keep);
yce  = yce(keep);

A    = 0.5 * abs(X(ncon(:,1)).*(Y(ncon(:,2))-Y(ncon(:,3))) + ...
                 X(ncon(:,2)).*(Y(ncon(:,3))-Y(ncon(:,1))) + ...
                 X(ncon(:,3)).*(Y(ncon(:,1))-Y(ncon(:,2))));
good = A > 1e-14 * sc^2;
ncon = ncon(good,:);
xce  = xce(good);
yce  = yce(good);

%% 8. REINDEX
used      = unique(ncon(:));
map       = zeros(numel(X),1);
map(used) = 1:numel(used);
X    = X(used);
Y    = Y(used);
ncon = map(ncon);

%% 9. MATERIAL  (1=body, 2=insert)
material = ones(size(ncon,1),1);
material(inpolygon(xce, yce, Xt, Yt)) = 2;

%% 10. REPORT
n_nodes    = numel(X);
n_elements = size(ncon,1);
fprintf('GenerateBevelToolMesh: %d nodes, %d elements\n', n_nodes, n_elements);
fprintf('  Material 1 (body): %d   Material 2 (insert): %d\n', ...
        sum(material==1), sum(material==2));
end
