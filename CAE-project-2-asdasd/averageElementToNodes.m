function nodalValues = averageElementToNodes(ncon, n_nodes, elemValues)

nodalValues = zeros(n_nodes,1);
count = zeros(n_nodes,1);

for i = 1:size(ncon,1)

    nodes = ncon(i,:);

    for j = 1:3
        n = nodes(j);

        nodalValues(n) = nodalValues(n) + elemValues(i);
        count(n) = count(n) + 1;
    end
end

% Avoid divide-by-zero if any isolated node exists
count(count == 0) = 1;

nodalValues = nodalValues ./ count;

end