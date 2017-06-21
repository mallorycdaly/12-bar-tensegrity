function q = findState_r(r, rod_pairs)
% This function finds the state vector of each rod given the nodal
% positions. Based on the coordinate system in Fig. 4 of Kyunam Kim's
% "Robust Learning of Tensegrity Robot Control for Locomotion through
% Form-Finding."
%
% The inputs are the following:
%   r: matrix of x,y,z, position of nodes
%   rod_pairs: each row defines the node indices corresponding to that rod
%
% The outputs are the following:
%   q: each row represents a rod, and the columns are the states
%       x,y,z,phi,theta, where x,y,z are the coordinates of the lower
%       (by z) node on the rod
%
% This should be general to any tensegrity, but I haven't tested on
% anything but a three-strut.

% Grab number of nodes and rods
num_nodes = size(r,1);
num_rods = num_nodes/2;

% Reorder nodes by rod pairs and so that higher (z) point is second
% Odd nodes: lower value of z
% Even nodes: higher value of z
r_rod = zeros(size(r));
for i = 1:num_rods
    if (r(rod_pairs(i,1),3) < r(rod_pairs(i,2),3))
        r_rod(2*i-1:2*i,:) = [r(rod_pairs(i,1),:)
            r(rod_pairs(i,2),:)];
    else
        r_rod(2*i-1:2*i,:) = [r(rod_pairs(i,2),:)
            r(rod_pairs(i,1),:)];
    end
end

% Find length of rod
rod_length = norm(r_rod(1,:)-r_rod(2,:));

% Find angles
phi = zeros(num_rods,1);
theta = zeros(num_rods,1);
for i = 1:num_rods
    % Displacement of x, y, and z of top node from bottom of rod
    x = r_rod(2*i,1) - r_rod(2*i-1,1);
    y = r_rod(2*i,2) - r_rod(2*i-1,2);
    z = r_rod(2*i,3) - r_rod(2*i-1,3);
    % Since we've already ordered the rods so that the greater value of z
    % is the first node, we'll have no trigonometric problems later on with
    % phi, which is measured from the z axis
    phi(i) = acosd(z/rod_length);
    % We will have trigonometric problems with finding theta. So we'll use
    % atan2 to solve them for us.
    theta(i) = atan2d(y,x);
end
q = [r_rod(1:2:num_nodes,:) phi theta];

% % Refind nodes based on q0 to check
% nodes_refind = zeros(size(nodes_rod));
% for i = 1:num_rods
%     nodes_refind(2*i-1,:) = nodes_rod(2*i-1,:);
%     nodes_refind(2*i,:) = nodes_rod(2*i-1,:) + ...
%         L*[sind(phi(i))*cosd(theta(i)) ...
%         sind(phi(i))*sind(theta(i)) ...
%         cosd(phi(i))];
% end
% nodes_rod
% nodes_refind