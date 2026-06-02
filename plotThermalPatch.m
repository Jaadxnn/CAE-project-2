function plotThermalPatch(ax, X, Y, ncon, T)

% =========================================================
% PLOTTHERMALPATCH
%
% Plots nodal temperature field using interpolated patch.
% =========================================================

cla(ax,'reset')
hold(ax,'on')

patch(ax, ...
    'Faces', ncon, ...
    'Vertices', [X Y], ...
    'FaceVertexCData', T, ...
    'FaceColor', 'interp', ...
    'EdgeColor', 'none');

axis(ax,'equal')
grid(ax,'on')

xlim(ax, [min(X) max(X)])
ylim(ax, [min(Y) max(Y)])

colormap(ax,'jet')

cb = colorbar(ax);
cb.Label.String = 'Temperature (°C)';

title(ax, 'Thermal Analysis: Temperature Distribution')
xlabel(ax, 'X')
ylabel(ax, 'Y')

hold(ax,'off')

end