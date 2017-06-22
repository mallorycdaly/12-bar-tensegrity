function F_cable = findCableForce(r, cable_pair, num_nodes, k_spring, ...
    L0_spring, rest_length)
% This function calculates the cable forces on each node. The cable is
% modeled as an inflexible string in series with a spring. The length of
% the string is the rest length. The force magnitude of the spring is found
% as F = k*(L-L0).
%
% The inputs are the following:
%   r: matrix of x,y,z, position of nodes
%   cable_pair: each row defines the node indices corresponding to that
%     cable
%   num_nodes: number of nodes
%   k_spring: spring constant of the springs
%   L0_spring: vector of initial lengths of the springs
%   rest_length: lengths of strings in the cable, defined as the cable rest
%       lengths, such that spring length can be found as total cable length
%       (node-to-node distance) minus rest length
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
    [cable,node_i_idx] = find(cable_pair == node_i);
    
    % Loop through each cable
    node_j_idx = zeros(size(node_i_idx));
    for n = 1:length(cable)
    
        % Find corresponding node for the cable
        if node_i_idx(n) == 1
            node_j_idx(n) = 2;
        else
            node_j_idx(n) = 1;
        end
        node_j = cable_pair(cable(n),node_j_idx(n));
        
        % Find cable length and direction vector
        cable_length = norm(r(node_j,:) - r(node_i,:));
        e_C = (r(node_j,:) - r(node_i,:)) / cable_length;
        
        % Cables can't push, so only add force is spring length is positive
        L_spring = cable_length - rest_length(cable(n));
        if L_spring > 0
            F_cable(i,:) = F_cable(i,:) + k_spring* ...
                (L_spring - L0_spring(cable(n))) * e_C;
        end
        
    end
end