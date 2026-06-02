function T = SolveThermal( ...
    n_nodes, n_elements, ...
    ncon, X, Y, ...
    t, ...
    material, ...
    k_tool, k_BUE, ...
    fixedNodes, fixedValues, ...
    convEdges, h, T_ambient)

if nargin < 8 || isempty(material)
    material = ones(n_elements,1);
end

if nargin < 13 || isempty(convEdges);  convEdges = [];  end
if nargin < 14 || isempty(h);          h = 0;           end
if nargin < 15 || isempty(T_ambient);  T_ambient = 0;   end

Kt = zeros(n_nodes, n_nodes);
Q  = zeros(n_nodes, 1);

%Assembly of conduction in the tool
for e = 1:n_elements

    n1 = ncon(e,1); n2 = ncon(e,2); n3 = ncon(e,3);
    x1 = X(n1); y1 = Y(n1);
    x2 = X(n2); y2 = Y(n2);
    x3 = X(n3); y3 = Y(n3);

    Ae = abs(0.5 * det([1 x1 y1; 1 x2 y2; 1 x3 y3]));

    if Ae <= 0
        error('Element %d has zero or negative area.', e);
    end

    b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
    c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;

    Bth = (1/(2*Ae)) * [b1 b2 b3; c1 c2 c3];

    if material(e) == 1
        k = k_tool;
    elseif material(e) == 2
        k = k_BUE;
    else
        error('Unknown thermal material ID in element %d', e);
    end

    Kte = t * Ae * Bth.' * (k * eye(2)) * Bth;

    nodes = [n1 n2 n3];
    for i = 1:3
        for j = 1:3
            Kt(nodes(i), nodes(j)) = Kt(nodes(i), nodes(j)) + Kte(i,j);
        end
    end

end

%The assembly of convection by extracting nodes from the edges found
if ~isempty(convEdges) && h > 0

    for ec = 1:size(convEdges, 1)

        ni = convEdges(ec,1);
        nj = convEdges(ec,2);
        L  = sqrt((X(nj)-X(ni))^2 + (Y(nj)-Y(ni))^2);

        Kce = (h*L*t/6)*[2 1; 1 2];
        Qce = (h*L*t*T_ambient/2)*[1; 1];

        nodes = [ni nj];
        for i = 1:2
            Q(nodes(i)) = Q(nodes(i)) + Qce(i);
            for j = 1:2
                Kt(nodes(i), nodes(j)) = Kt(nodes(i), nodes(j)) + Kce(i,j);
            end
        end

    end

end


% DIRICHLET BCS or fixed temepratures application
fixedNodes  = fixedNodes(:);
fixedValues = fixedValues(:);

if length(fixedNodes) ~= length(fixedValues)
    error('fixedNodes and fixedValues must have the same length.');
end

for i = 1:length(fixedNodes)

    node = fixedNodes(i);
    Tval = fixedValues(i);

    Q = Q - Kt(:,node) * Tval;
    Kt(node,:) = 0;  Kt(:,node) = 0;
    Kt(node,node) = 1;
    Q(node) = Tval;

end

if rcond(Kt) < 1e-14
    warning('Thermal matrix is close to singular. Check thermal boundary conditions.');
end

T = Kt \ Q;

end