function [F_rod, L_rod] = findRodForce(r, rod_pair, num_nodes, ...
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

for i = 1:num_rods
    
    % Define node i as lower node and k as higher node
    if r(rod_pair(i,1),3) < r(rod_pair(i,2),3)
        node_i = rod_pair(i,1);
        node_k = rod_pair(i,2);
    else
        node_i = rod_pair(i,2);
        node_k = rod_pair(i,1);
    end
    
    % Calculate rod length and direction along rod
    L_rod(i) = norm(r(node_k,:) - r(node_i,:));     
    e_R = (r(node_k,:) - r(node_i,:)) / L_rod(i);
    
    % Calculate rod force magnitude
    F_ik = k_rod*(L_rod(i) - L0_rod);
    
    % Assign force values
    F_rod(node_i,:) = F_ik*e_R;
    F_rod(node_k,:) = -F_ik*e_R;
    
end