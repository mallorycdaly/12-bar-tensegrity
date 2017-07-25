function [F_rod, L_rod] = findRodForce2(r, rod_pair, num_nodes, ...
    num_rods, k_rod, L0_rod)
% This function calculates the rod forces on each node. The rod can be
% thought of as a very stiff spring. The force magnitude is calculated as
% F = k*(L-L0).
% 
% The inputs are the following: 
%   r: matrix of x,y,z, position of nodes
%   rod_pair: each row defines the node indices corresponding to that rod
%   num_nodes: number of nodes
%   num_rods: number of rods
%   k_rod: spring constant of the rod (could be found as EA/L)
%   L0_rod: length of rods before deformation
%
% The outputs are the following:
%   F_rod: each row defines the force vector (x,y,z) acting on the node
%       from the rod
%   L_rod: each row defines the rod's length

% Initialize variables
F_rod = zeros(num_nodes,3);
L_rod = zeros(num_rods,1);

% Find sum of forces on node from each cable
for i = 1:num_nodes
    
    % Find node connected to current node by a rod
    node_i = i;
    [cable,node_i_idx] = find(rod_pair == node_i);

    % Find corresponding node for the rod
    if node_i_idx == 1
        node_k_idx = 2;
    else
        node_k_idx = 1;
    end
    node_k = rod_pair(cable,node_k_idx);
        
    % Find rod length and direction vector
    L_rod(i) = norm(r(node_k,:) - r(node_i,:));
    e_R = (r(node_k,:) - r(node_i,:)) / L_rod(i);

    % Calculate rod force
    F_rod(i,:) = k_rod*(L_rod(i) - L0_rod) * e_R;
        
end