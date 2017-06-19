function [a, F_rod, F] = findForcesAndAccelerations(r, q, q_dot, ...
    F_cable, rod_pairs, m)
% This function calculates the cable forces on each node.
%
% The inputs are the following:
%   r: matrix of x,y,z, position of nodes
%   v: matrix of x,y,z, velocity of nodes
%   cable_pairs: each row defines the node indices corresponding to that
%     cable
%   spring_constant: spring constant
%   rest_lengths: rest length of cable, such that spring length can be
%       found as total cable length (node-to-node distance) - rest length
%
% The outputs are the following:
%   F_cable: each row defines the force vector (x,y,z) acting on the node
%       from the cables

% Grab number of nodes
num_nodes = size(r,1);
num_rods = size(rod_pairs,1);

% Find accelerations and forces
a = zeros(num_nodes,3);
F_rod = zeros(num_nodes,3);
F = zeros(num_nodes,3);
for i = 1:num_rods
    
    % To do: Could pass this information in from previous functions, or
    % make a function that does this
    % Define node i as lower node and k as higher node
    if r(rod_pairs(i,1),3) < r(rod_pairs(i,2),3)
        node_i = rod_pairs(i,1);
        node_k = rod_pairs(i,2);
    else
        node_i = rod_pairs(i,2);
        node_k = rod_pairs(i,1);
    end
    
    % Find rod direction
    L = norm(r(node_k,:) - r(node_i,:));  % To do: Check if L is constant  
    e_R = ((r(node_k,:) - r(node_i,:)) / L)';
    
    % Extract needed values
    phi = q(i,4);
    theta = q(i,5);
    phi_dot = q_dot(i,4);
    theta_dot = q_dot(i,5);
    
    % Find accelerations knowing that m*a = F_cable(r) + F_rod(r,v,a)
    % To do: What happens if M is noninvertible?
    
        % Set up and solve M*a = b for i nodes
        C1 = sind(phi)*cosd(theta);
        C2 = sind(phi)*sind(theta);
        C3 = cosd(phi);
        M = [(1+C1*e_R(1))     C2*e_R(1)     C3*e_R(1);
                 C1*e_R(2) (1+C2*e_R(2))     C3*e_R(2);
                 C1*e_R(3)     C2*e_R(3) (1+C3*e_R(1))];
        b = 1/m*( (m*(L*phi_dot^2 + L*theta_dot^2*sind(phi)^2) + ...
            F_cable(node_k,:)*e_R)*e_R + F_cable(node_i,:)');
        a(node_i,:) = (M\b)';
        
        % Set up and solve M*a = b for k nodes
        C1 = sind(phi)*cosd(theta);
        C2 = sind(phi)*sind(theta);
        C3 = cosd(phi);
        M = [(1-C1*e_R(1))    -C2*e_R(1)    -C3*e_R(1);
                -C1*e_R(2) (1-C2*e_R(2))    -C3*e_R(2);
                -C1*e_R(3)    -C2*e_R(3) (1-C3*e_R(1))];
        b = 1/m*( -(m*(L*phi_dot^2 + L*theta_dot^2*sind(phi)^2) + ...
            F_cable(node_k,:)*e_R)*e_R + F_cable(node_k,:)');
        a(node_k,:) = (M\b)';
        
        % Check by finding forces
%         F_rod(node_i,:) = m*(a(node_i,1)*C1 + a(node_i,2)*C2 + ...
%             a(node_i,3)*C3 - L*phi_dot

end
