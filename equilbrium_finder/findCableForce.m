function F_cable = findCableForce(r, cable_pair, num_nodes, ...
    spring_constant, spring_initial_length, rest_length)
% This function calculates the cable forces on each node. The cable is
% modeled as an inflexible string in series with a spring. The length of
% the string is the rest length. The force in the spring is found as
%   F = k*(l-l_0).
%
% The inputs are the following:
%   r: matrix of x,y,z, position of nodes
%   cable_pair: each row defines the node indices corresponding to that
%     cable
%   num_nodes: number of nodes
%   spring_constant: spring constant
%   spring_initial_length: initial length (l_0) of the spring
%   rest_length: rest lengths of cable, such that spring length can be
%       found as total cable length (node-to-node distance) minus rest
%       length
%
% The outputs are the following:
%   F_cable: each row defines the force vector (x,y,z) acting on the node
%       from the cables

% Initialize cable force matrix
F_cable = zeros(num_nodes,3);

% Find sum of forces on node from each cable
for i = 1:num_nodes
    
    % Find nodes connected to current node by cables
    node_i = i;
    [row_ij,col_i] = find(cable_pair == node_i);
    
    % Loop through each cable
    col_j = zeros(size(col_i));
    for n = 1:size(row_ij)
    
        % Find corresponding node for the cable
        if col_i(n) == 1
            col_j(n) = 2;
        else
            col_j(n) = 1;
        end
        node_j = cable_pair(row_ij(n),col_j(n));
        
        % Calculate cable force and add to total cable force for this node
        cable_length = norm(r(node_j,:) - r(node_i,:));
        e_C = (r(node_j,:) - r(node_i,:)) / cable_length;
        % Cables can't push, so only add force is spring length is positive
        spring_length = cable_length - rest_length(i);
        if spring_length > 0
            F_cable(i,:) = F_cable(i,:) + spring_constant* ...
                (spring_length - spring_initial_length) * e_C;
        end
        
    end
end