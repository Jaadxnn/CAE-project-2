function T = SolveThermalGeneral( ...
    n_nodes, n_elements, ...
    ncon, X, Y, ...
    t, ...
    kx, ky, ...
    Q_uniform, ...
    fixedNodes, fixedValues, ...
    convectionNodes, h, T_ambient, nx, ny)

% Default optional args
if nargin < 13 || isempty(convectionNodes); convectionNodes = []; end
if nargin < 14 || isempty(h);               h = 0;                end
if nargin < 15 || isempty(T_ambient);       T_ambient = 0;        end

Kt = zeros(n_nodes, n_nodes);
Q  = zeros(n_nodes, 1);

%Assemble gobal stiffness and heat source matrices
for e = 1:n_elements
    n1 = ncon(e,1);
    n2 = ncon(e,2);
    n3 = ncon(e,3);

    x1 = X(n1); y1 = Y(n1);
    x2 = X(n2); y2 = Y(n2);
    x3 = X(n3); y3 = Y(n3);

    Ae = 0.5 * abs(det([1 x1 y1;
                         1 x2 y2;
                         1 x3 y3]));

    if Ae <= 0
        error('Element %d has zero or negative area.', e);
    end

    b1 = y2 - y3;  b2 = y3 - y1;  b3 = y1 - y2;
    c1 = x3 - x2;  c2 = x1 - x3;  c3 = x2 - x1;

    %Thermal matrix
    Bth = (1/(2*Ae)) * [b1 b2 b3;
                         c1 c2 c3];

    %Anisotropic conductivity matrix
    Dth = [kx 0;
            0 ky];

    %Element stiffness matrix
    Kte = t * Ae * (Bth.') * Dth * Bth;

    %Uniform heat source contribution that is distributed evenly to 3 nodes
    Qe = (Q_uniform * Ae * t / 3) * ones(3,1);

    %Assemble
    nodes = [n1 n2 n3];
    for i = 1:3
        Q(nodes(i)) = Q(nodes(i)) + Qe(i);
        for j = 1:3
            Kt(nodes(i), nodes(j)) = Kt(nodes(i), nodes(j)) + Kte(i,j);
        end
    end
end


%The assembly of convection BC using the Robin method
if ~isempty(convectionNodes) && h > 0
    convectionNodes = sort(convectionNodes(:));
    
    for i = 1:length(convectionNodes)-1
        ni = convectionNodes(i);
        nj = convectionNodes(i+1);

        xi = X(ni); yi = Y(ni);
        xj = X(nj); yj = Y(nj);
        L  = sqrt((xj-xi)^2 + (yj-yi)^2);

        %Skip if the two nodes are not actually adjacent to each other
        %In this case if their distance is more than 1.5x the expected edge length
        dx_expected = abs(max(X) - min(X)) / (nx - 1);
        dy_expected = abs(max(Y) - min(Y)) / (ny - 1);
        L_max = 1.5 * max(dx_expected, dy_expected);

        if L > L_max
            continue
        end

        Kce = (h * L * t / 6) * [2 1; 1 2];
        Qce = (h * L * t * T_ambient / 2) * [1; 1];

        idx = [ni nj];
        for a = 1:2
            Q(idx(a)) = Q(idx(a)) + Qce(a);
            for b = 1:2
                Kt(idx(a), idx(b)) = Kt(idx(a), idx(b)) + Kce(a,b);
            end
        end
    end
end


%Apply the fixed temperature BC (dirichelet)
fixedNodes  = fixedNodes(:);
fixedValues = fixedValues(:);

if length(fixedNodes) ~= length(fixedValues)
    error('fixedNodes and fixedValues must have the same length.');
end

for i = 1:length(fixedNodes)
    node = fixedNodes(i);
    Tval = fixedValues(i);

    %Modify the RHS for non zero Dirichlet
    Q = Q - Kt(:, node) * Tval;

    %Zero out the rows and columns, and set the diagonal to 1
    Kt(node, :) = 0;
    Kt(:, node) = 0;
    Kt(node, node) = 1;
    Q(node) = Tval;
end
%Solve and checking if the matrix has no tbeen formed properly
if rcond(Kt) < 1e-14
    warning('Thermal stiffness matrix is near singular. Check BCs.');
end

T = Kt \ Q;
end