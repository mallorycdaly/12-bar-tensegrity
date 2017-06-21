function [r, cable_pair, rod_pair, num_nodes, num_cables, num_rods, ...
    L_cable, L_rod] = formThreeBar(edge_length, height, rotation_angle)
% This function forms a three strut tensegrity prism using the numbering
% scheme in Fig. 10. of "Review of Form-Finding Methods for Tensegrity
% Structures" by A. Tibert.
%
% The inputs are the following:
%   edge_length: the edge length of the top and bottom triangles
%   height: the distance between the top and bottom triangles
%   rotation_angle: the rotation angle between the top and bottom triangles

%
% The outputs are the following:
%   r: matrix of x,y,z, position of nodes based on design parameters
%   cable_pair: each row defines the node indices corresponding to that
%     cable
%   rod_pair: each row defines the node indices corresponding to that rod
%   num_nodes: number of nodes
%   num_cables: number of cables
%   num_rods: number of rods
%   L_cable: lengths of the cables (node-to-node distance)
%   L_rod: lengths of the rods

% Base triangle centered at (0,0,0)
node4 = [edge_length/2 -sqrt(3)/6*edge_length 0];
node5 = [0 sqrt(3)/3*edge_length 0];
node6 = [-edge_length/2 -sqrt(3)/6*edge_length 0];

% Rotation matrix
rotM = [cosd(rotation_angle) -sind(rotation_angle) 0;
        sind(rotation_angle)  cosd(rotation_angle) 0;
                           0                     0 1];

% Rotate base triangle and translate to form top triangle
node1 = node4*rotM + [0 0 height];
node2 = node5*rotM + [0 0 height];
node3 = node6*rotM + [0 0 height];

% Assemble nodes in matrix
r = [node1; node2; node3; node4; node5; node6];

% Cables by pairs of node indices
cable_pair = [1 2;     % cable 1
               2 3;     % cable 2
               1 3;     % cable 3
               1 4;     % cable 4
               2 5;     % cable 5
               3 6;     % cable 6
               4 5;     % cable 7
               5 6;     % cable 8
               4 6];    % cable 9
           
% Rods by pairs of node indices           
rod_pair   = [1 6;     % bar 10
               2 4;     % bar 11
               3 5];    % bar 12

% Find number of nodes, cables, and rods           
num_nodes = size(r,1);
num_cables = size(cable_pair,1);
num_rods = size(rod_pair,1);

% Original cable and rod lengths, before any deformation
L_cable = zeros(num_cables,1);
for i = 1: num_cables
    L_cable(i) = norm(r(cable_pair(i,1),:) - ...
        r(cable_pair(i,2),:));
end
L_rod = norm(r(rod_pair(1,1),:) - r(rod_pair(1,2),:));