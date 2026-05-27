function [Xdef, Ydef, Xorig, Yorig] = ...
    display_structure(n_element, ncon, X, Y, U, scale)

Xdef  = zeros(n_element,4);
Ydef  = zeros(n_element,4);

Xorig = zeros(n_element,4);
Yorig = zeros(n_element,4);

for i = 1:n_element

    n1 = ncon(i,1);
    n2 = ncon(i,2);
    n3 = ncon(i,3);

    % Original coordinates
    x1 = X(n1); y1 = Y(n1);
    x2 = X(n2); y2 = Y(n2);
    x3 = X(n3); y3 = Y(n3);

    % Store original closed triangle
    Xorig(i,:) = [x1 x2 x3 x1];
    Yorig(i,:) = [y1 y2 y3 y1];

    % Displacements
    u1 = U(2*n1-1);
    v1 = U(2*n1);

    u2 = U(2*n2-1);
    v2 = U(2*n2);

    u3 = U(2*n3-1);
    v3 = U(2*n3);

    % Deformed coordinates
    Xdef(i,:) = [ ...
        x1 + scale*u1,...
        x2 + scale*u2,...
        x3 + scale*u3,...
        x1 + scale*u1];

    Ydef(i,:) = [ ...
        y1 + scale*v1,...
        y2 + scale*v2,...
        y3 + scale*v3,...
        y1 + scale*v1];

end

end