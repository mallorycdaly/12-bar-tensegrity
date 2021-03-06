function [q, q_dot] = findState_rv(r, v, rod_pairs)
% This function finds the derivative of the state vector of each rod given
% the nodal positions. Based on the coordinate system in Fig. 4 of Kyunam
% Kim's "Robust Learning of Tensegrity Robot Control for Locomotion through
% Form-Finding."
%
% The inputs are the following:
%   r: matrix of x,y,z, position of nodes
%   v: matrix of x,y,z velocities of nodes
%   rod_pairs: each row defines the node indices corresponding to that rod
%
% The outputs are the following:
%   q_dot: each row represents a rod, and the columns are the states
%       x_dot,y_dot,z_dot,phi_dot,theta_dot, in the same rod order as q
%
% This should be general to any tensegrity, but I haven't tested on
% anything but a three-strut.

%% Find states

% Grab number of nodes and rods
num_nodes = size(r,1);
num_rods = num_nodes/2;

% Reorder r by rod pairs and so that higher (z) point is second and store
% vector of ordered node indices
% Odd nodes: lower value of z
% Even nodes: higher value of z
rod_pairs_ordered = zeros(num_nodes,1);
r_rod = zeros(size(r));
for i = 1:num_rods
    if (r(rod_pairs(i,1),3) < r(rod_pairs(i,2),3))
        r_rod(2*i-1:2*i,:) = [r(rod_pairs(i,1),:)
            r(rod_pairs(i,2),:)];
        rod_pairs_ordered(2*i-1:2*i) = [rod_pairs(i,1); rod_pairs(i,2)];
    else
        r_rod(2*i-1:2*i,:) = [r(rod_pairs(i,2),:)
            r(rod_pairs(i,1),:)];
        rod_pairs_ordered(2*i-1:2*i) = [rod_pairs(i,2); rod_pairs(i,1)];
    end
end

% Reorder v to match order of r_rod
v_rod = zeros(size(r));
for i = 1:num_nodes
    v_rod(i,:) = v(rod_pairs_ordered(i),:);
end

% Find length of rod
L = norm(r_rod(1,:)-r_rod(2,:));

% Find angles
phi = zeros(num_rods,1);
theta = zeros(num_rods,1);
phi_dot = zeros(num_rods,1);
theta_dot = zeros(num_rods,1);
for i = 1:num_rods
    
    % Displacement of x, y, and z of top node from bottom of rod
    del_x = r_rod(2*i,1) - r_rod(2*i-1,1);
    del_y = r_rod(2*i,2) - r_rod(2*i-1,2);
    del_z = r_rod(2*i,3) - r_rod(2*i-1,3);
    
    % Since we've already ordered the rods so that the greater value of z
    % is the first node, we'll have no trigonometric problems later on with
    % phi, which is measured from the z axis
    phi(i) = acosd(del_z/L);
    
    % We will have trigonometric problems with finding theta. So we'll use
    % atan2 to solve them for us.
    theta(i) = atan2d(del_y,del_x);
    
    % Extract x_dot and z_dot at each node to make things easier
    x_dot_i = v_rod(2*i-1,1);
    x_dot_k = v_rod(2*i,1);
    z_dot_i = v_rod(2*i-1,3);
    z_dot_k = v_rod(2*i,3);
    
    % Find phi_dot, and check for singularity
    if sind(phi(i)) == 0
        phi_dot(i) = 0;
    else
        phi_dot(i) = 1/(-L*sind(phi(i)))*(z_dot_k - z_dot_i);
    end
    
    % Find theta_dot, and check for singularity
    if sind(phi(i))*sind(theta(i)) == 0
        theta_dot(i) = 0;
    else
        theta_dot(i) = 1/(-sind(phi(i))*sind(theta(i))) * ...
            (1/L*(x_dot_k - x_dot_i) - phi_dot(i)*cosd(phi(i))* ...
            cosd(theta(i)));
    end
    
end
q = [r_rod(1:2:num_nodes,:) phi theta];
q_dot = [v_rod(1:2:num_nodes,:) phi_dot theta_dot];

%% Check for errors
error_tol = 1e-6;

% Refind r_rod based on q to check
r_rod_refind = zeros(size(r_rod));
for i = 1:num_rods
    r_rod_refind(2*i-1,:) = r_rod(2*i-1,:);
    r_rod_refind(2*i,:) = r_rod(2*i-1,:) + ...
        L*[sind(phi(i))*cosd(theta(i)) ...
           sind(phi(i))*sind(theta(i)) ...
           cosd(phi(i))];
end

% Compare and throw error if there is mismatch
if (any(any(r_rod - r_rod_refind > error_tol)))
    error(['There is a mismatch between the nodal positions and the ' ...
        'reconstructed nodal positions based on the calculated angles.'])
end

% Refind v_rod based on q_dot to check
v_rod_refind = zeros(size(r_rod));
for i = 1:num_rods
    v_rod_refind(2*i-1,:) = v_rod(2*i-1,:);
    v_rod_refind(2*i,:) = v_rod(2*i-1,:) + ...
        L*[phi_dot(i)*cosd(phi(i))*cosd(theta(i)) - ...
               theta_dot(i)*sind(phi(i))*sind(theta(i)) ...
           phi_dot(i)*cosd(phi(i))*sind(theta(i)) - ...
               theta_dot(i)*sind(phi(i))*cosd(theta(i)) ...
           -phi_dot(i)*sind(phi(i))];
end

% Compare and throw error if there is mismatch
if (any(any(v_rod - v_rod_refind > error_tol)))
    error(['There is a mismatch between the nodal velocities and the ' ...
        'reconstructed nodal velocities based on the calculated angles.'])
end
