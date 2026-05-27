clc;
clear;
close all;

%% Read Excel File
data = readmatrix('newansysmesh_converted.xlsx');

%% Connectivity Matrix
% Columns:
% ncon1 = col 3
% ncon2 = col 4
% ncon3 = col 5

ncon = data(:,3:5);

% Remove NaN rows
ncon = ncon(~any(isnan(ncon),2),:);

%% Node Coordinates
% X = col 6
% Y = col 7

coords = data(:,6:7);

% Remove NaN rows
coords = coords(~any(isnan(coords),2),:);

X = coords(:,1);
Y = coords(:,2);

%% Plot Mesh
figure;
hold on;
axis equal;
box on;

for e = 1:size(ncon,1)

    % Node numbers for current element
    nodes = ncon(e,:);

    % Element coordinates
    xe = X(nodes);
    ye = Y(nodes);

    % Close triangle
    xe = [xe; xe(1)];
    ye = [ye; ye(1)];

    % Plot element
    plot(xe, ye, 'k-');

end

%% Plot nodes
scatter(X, Y, 15, 'filled');

xlabel('X Coordinate');
ylabel('Y Coordinate');
title('ANSYS Mesh');

legend('Elements','Nodes');