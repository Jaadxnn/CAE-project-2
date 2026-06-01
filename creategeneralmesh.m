function [X,Y,ncon, n_nodes, n_elements] = creategeneralmesh(X_g,Y_g, nx, ny)

%Top edge
xL = X_g(1);
xR = X_g(4);
yL = Y_g(1);
yR = Y_g(4);
x_top = linspace(xL,xR,nx);
y_top = linspace(yL,yR, nx);

%Bottom edge
xL_bot = X_g(2);
xR_bot = X_g(3);
yL_bot = Y_g(2);
yR_bot = Y_g(3);
x_bot = linspace(xL_bot,xR_bot, nx);
y_bot = linspace(yL_bot,yR_bot, nx);



nx = length(x_top);
t_vals = linspace(0,1,ny);
X = [];
Y = [];

% Transfinite interpolation
for j = 1:ny
    tau = t_vals(j);
    x_row = (1 - tau) * x_top + tau * x_bot;
    y_row = (1 - tau) * y_top + tau * y_bot;
    idx = (j-1)*nx + (1:nx);
    X(idx) = x_row;
    Y(idx) = y_row;
end

% Connectivity with CCW winding
ncon = zeros(2*(nx-1)*(ny-1), 3);
e = 0;
for j = 1:(ny-1)
    for i = 1:(nx-1)
        n1 = (j-1)*nx + i;
        n2 = n1 + 1;
        n3 = n1 + nx;
        n4 = n3 + 1;
        e = e + 1; ncon(e,:) = [n1 n3 n2];  % CCW
        e = e + 1; ncon(e,:) = [n2 n3 n4];  % CCW
    end
end

n_nodes    = length(X);
n_elements = size(ncon, 1);


end