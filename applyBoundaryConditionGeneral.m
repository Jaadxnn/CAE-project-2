function [fixedNodes, fixedValues, insulatedNodes, convectionNodes] = ...
    applyBoundaryConditionGeneral(side, bcType, nx, ny, fixedTemp)

    % Get nodes for the requested side
    switch side
        case 'bottom'
            nodes = (ny-1)*nx + (1:nx);
        case 'top'
            nodes = 1:nx;
        case 'left'
            nodes = 1:nx:(ny-1)*nx + 1;
        case 'right'
            nodes = nx:nx:ny*nx;
    end

    % Initialise outputs
    fixedNodes      = [];
    fixedValues     = [];
    insulatedNodes  = [];
    convectionNodes = [];

    % Apply the BC type
    switch bcType
        case 'Fixed Temperature'
            fixedNodes  = nodes(:);
            fixedValues = ones(length(nodes), 1) * fixedTemp;

        case 'Insulated'
            insulatedNodes = nodes(:);

        case 'Convection'
            convectionNodes = nodes(:);
    end
end