function F_rod = findRodForce(r, rod_pairs, q, q_dot, q_dotdot, ...
    F_cable, mass)
% This function calculates the rod forces on each node.
% 
% The inputs are the following:
%   r: matrix of x,y,z, position of nodes
%   rod_pairs: each row defines the node indices corresponding to that rod
%   q: state vector x,y,z,phi,theta
%   q_dot: first derivative of state vector with respect to time
%   q_dotdot: second derivative of state vector with respect to time
%   F_cable: each row defines the force vector (x,y,z) acting on the node
%       from the cables
%   mass: the mass of each node
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

    % Extract needed values
    x_dotdot = q_dotdot(i,1);
    y_dotdot = q_dotdot(i,2);
    z_dotdot = q_dotdot(i,3);
    phi = q(i,4);
    theta = q(i,5);
    phi_dot = q_dot(i,4);
    theta_dot = q_dot(i,5);

    % Find rod's force magnitude
    rod_length = norm(r(node_k,:) - r(node_i,:));     
    e_R = (r(node_k,:) - r(node_i,:)) / rod_length;
    F_ik = mass*(x_dotdot*sind(phi)*cosd(theta) + ...
        y_dotdot*sind(phi)*sind(theta) + z_dotdot*cosd(phi) - ...
        rod_length*phi_dot^2 - rod_length*theta_dot^2*sind(phi)^2) - ...
        F_cable(node_k,:)*e_R';
    
    % Assign force to each node of rod
    F_rod(node_i,:) = -F_ik*e_R;
    F_rod(node_k,:) = F_ik*e_R;
    
end