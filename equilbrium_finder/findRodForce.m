function F_rod = findRodForce(r, rod_pairs, spring_constant, ...
    original_rod_length)
% This function calculates the rod forces on each node. The rod can be
% thought of as a very stiff spring. The force is calculated as
%   F = -spring_constant * (rod_length - original_rod_length)
% 
% The inputs are the following:
%   r: matrix of x,y,z, position of nodes
%   rod_pairs: each row defines the node indices corresponding to that rod
%   spring_constant: spring constant of the rod (could be found as EA/L)
%
% The outputs are the following:
%   F_rod: each row defines the force vector (x,y,z) acting on the node
%       from the rod

% Grab number of nodes and rods
num_nodes = size(r,1);
num_rods = size(rod_pairs,1);

% Initialize rod force matrix
F_rod = zeros(num_nodes,3);

for i = 1:num_rods
    
    % Define node i as lower node and k as higher node
    if r(rod_pairs(i,1),3) < r(rod_pairs(i,2),3)
        node_i = rod_pairs(i,1);
        node_k = rod_pairs(i,2);
    else
        node_i = rod_pairs(i,2);
        node_k = rod_pairs(i,1);
    end
    
    % Calculate rod length and direction along rod
    rod_length = norm(r(node_k,:) - r(node_i,:));     
    e_R = (r(node_k,:) - r(node_i,:)) / rod_length;
    
    % Calculate rod force magnitude
    F_ik = spring_constant*(rod_length - original_rod_length);
    
    % Assign force values
    F_rod(node_i,:) = -F_ik*e_R;
    F_rod(node_k,:) = F_ik*e_R;
    
end