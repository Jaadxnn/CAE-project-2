function plotSelectedContour(ax, dropdownValue, data)

% =========================================================
% PLOTSELECTEDCONTOUR
%
% Plots FEM contour lines based on dropdown selection.
%
% INPUTS:
%
% ax             -> UIAxes handle
% dropdownValue  -> selected dropdown string
% data           -> FEM data structure
%
% =========================================================

% =========================================================
% SELECT FIELD
% =========================================================
switch dropdownValue

    % -----------------------------------------------------
    % DISPLACEMENTS
    % -----------------------------------------------------
    case 'Ux'
        nodalField = data.Ux;

    case 'Uy'
        nodalField = data.Uy;

    case 'Uxy'
        nodalField = data.Umag;

    % -----------------------------------------------------
    % STRAINS
    % -----------------------------------------------------
    case 'Ex'
        nodalField = averageElementToNodes( ...
            data.ncon, ...
            data.n_nodes, ...
            data.Ex);

    case 'Ey'
        nodalField = averageElementToNodes( ...
            data.ncon, ...
            data.n_nodes, ...
            data.Ey);

    case 'Exy'
        nodalField = averageElementToNodes( ...
            data.ncon, ...
            data.n_nodes, ...
            data.Gxy);

    % -----------------------------------------------------
    % STRESSES
    % -----------------------------------------------------
    case 'Sx'
        nodalField = averageElementToNodes( ...
            data.ncon, ...
            data.n_nodes, ...
            data.Sx);

    case 'Sy'
        nodalField = averageElementToNodes( ...
            data.ncon, ...
            data.n_nodes, ...
            data.Sy);

    case 'Sxy'
        nodalField = averageElementToNodes( ...
            data.ncon, ...
            data.n_nodes, ...
            data.Sxy);

    otherwise
        error('Unknown contour selection');

end

% =========================================================
% DRAW CONTOURS
% =========================================================
drawContourLines( ...
    ax, ...
    data.X, ...
    data.Y, ...
    data.ncon, ...
    nodalField, ...
    dropdownValue);

end


% =========================================================
% ELEMENT -> NODE AVERAGING
% =========================================================
function nodalValues = averageElementToNodes( ...
    ncon, n_nodes, elemValues)

nodalValues = zeros(n_nodes,1);
count = zeros(n_nodes,1);

for i = 1:size(ncon,1)

    nodes = ncon(i,:);

    for j = 1:3

        n = nodes(j);

        nodalValues(n) = ...
            nodalValues(n) + elemValues(i);

        count(n) = count(n) + 1;

    end
end

nodalValues = nodalValues ./ count;

end


% =========================================================
% DRAW CONTOUR LINES
% =========================================================
function drawContourLines( ...
    ax, X, Y, ncon, nodalField, plotTitle)

cla(ax);
hold(ax,'on');

colormap(ax, jet);

cmap = colormap(ax);
ncol = size(cmap,1);

% =========================================================
% DRAW MESH
% =========================================================
patch(ax, ...
    'Faces', ncon, ...
    'Vertices', [X Y], ...
    'FaceColor', 'none', ...
    'EdgeColor', [0.7 0.7 0.7], ...
    'LineWidth', 0.5);

% =========================================================
% CONTOUR LEVELS
% =========================================================
n_levels = 10;

levels = linspace( ...
    min(nodalField), ...
    max(nodalField), ...
    n_levels);

% =========================================================
% MARCHING TRIANGLES
% =========================================================
for L = 1:length(levels)

    levelValue = levels(L);

    cidx = round( ...
        (L-1)/(length(levels)-1) * (ncol-1)) + 1;

    color = cmap(cidx,:);

    for i = 1:size(ncon,1)

        nodes = ncon(i,:);

        x = X(nodes);
        y = Y(nodes);

        s = nodalField(nodes);

        pts = [];

        edges = [1 2;
                 2 3;
                 3 1];

        for e = 1:3

            nA = edges(e,1);
            nB = edges(e,2);

            sA = s(nA);
            sB = s(nB);

            if (sA - levelValue) * ...
               (sB - levelValue) < 0

                lambda = ...
                    (levelValue - sA) / (sB - sA);

                x_cross = ...
                    x(nA) + lambda * (x(nB) - x(nA));

                y_cross = ...
                    y(nA) + lambda * (y(nB) - y(nA));

                pts = [pts;
                       x_cross y_cross];
            end
        end

        if size(pts,1) == 2

            plot(ax, ...
                pts(:,1), ...
                pts(:,2), ...
                'Color', color, ...
                'LineWidth', 1.5);

        end
    end
end

% =========================================================
% FORMATTING
% =========================================================
caxis(ax, ...
    [min(nodalField) max(nodalField)]);

cb = colorbar(ax);
cb.Label.String = plotTitle;

axis(ax,'equal');

title(ax, ['Contour Plot: ' plotTitle]);

xlabel(ax,'X');
ylabel(ax,'Y');

hold(ax,'off');

end