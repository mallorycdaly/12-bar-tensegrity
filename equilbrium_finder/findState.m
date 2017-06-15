function q = findState(nodes, bar_pairs)
% This function finds the state vector of each rod given the nodal
% positions. Based on the coordinate system in Fig. 4 of Kyunam Kim's
% "Robust Learning of Tensegrity Robot Control for Locomotion through
% Form-Finding."
%
% The inputs are the following:
%   nodes: matrix of x,y,z, position of nodes
%   bar_pairs: each row defines the node indices corresponding to that bar
%
% The outputs are the following:
%   q: each row represents a bar, and the columns are the states
%       x,y,z,phi,theta, where x,y,z are the coordinates of the lower
%       (by z) node on the rod
%
% This should be general to any tensegrity, but I haven't tested on
% anything but a three-strut.

% The nodes should be size numNodes x 3 (x,y,z)
num_nodes = size(nodes,1);
num_bars = num_nodes/2;

% Reorder nodes by bar pairs and so that higher (z) point is second
% Odd nodes: lower value of z
% Even nodes: higher value of z
nodes_bar = zeros(size(nodes));
for i = 1:num_bars
    if (nodes(bar_pairs(i,1),3) > nodes(bar_pairs(i,2),3))
        nodes_bar(2*i-1:2*i,:) = [nodes(bar_pairs(i,2),:)
            nodes(bar_pairs(i,1),:)];
    else
        nodes_bar(2*i-1:2*i,:) = [nodes(bar_pairs(i,1),:)
            nodes(bar_pairs(i,2),:)];
    end
end

% Find length of bar
L = norm(nodes_bar(1,:)-nodes_bar(2,:));

% Find angles
phi = zeros(num_bars,1);
theta = zeros(num_bars,1);
for i = 1:num_bars
    % Displacement of x, y, and z of top node from bottom of rod
    x = nodes_bar(2*i,1) - nodes_bar(2*i-1,1);
    y = nodes_bar(2*i,2) - nodes_bar(2*i-1,2);
    z = nodes_bar(2*i,3) - nodes_bar(2*i-1,3);
    % Since we've already ordered the rods so that the greater value of z
    % is the first node, we'll have no trigonometric problems later on with
    % phi, which is measured from the z axis
    phi(i) = acosd(z/L);
    % We will have trigonometric problems with finding theta. So we'll use
    % atan2 to solve them for us.
    theta(i) = atan2d(y,x);
end
q = [nodes_bar(1:2:num_nodes,:) phi theta];

% % Refind nodes based on q0 to check
% nodes_refind = zeros(size(nodes_bar));
% for i = 1:num_bars
%     nodes_refind(2*i-1,:) = nodes_bar(2*i-1,:);
%     nodes_refind(2*i,:) = nodes_bar(2*i-1,:) + ...
%         L*[sind(phi(i))*cosd(theta(i)) ...
%         sind(phi(i))*sind(theta(i)) ...
%         cosd(phi(i))];
% end
% nodes_bar
% nodes_refind