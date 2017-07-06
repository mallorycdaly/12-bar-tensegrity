function plotRotatedTensegrity(r, cable_pair, rod_pair, ground_face, ...
    num_cables_used, labels_on, cable_format)
% Plots tensegrity (likely one that has been actuated) so that its ground
% face is oriented downwards.
%
% The inputs are the following:
%   r: matrix of x,y,z, position of nodes
%   cable_pair: each row defines the node indices corresponding to that
%     cable
%   rod_pair: each row defines the node indices corresponding to that bar
%   ground_face: node indices of one ground face
%   num_cables_used: number of actuated cables
%   labels_on: boolean to add labels of bars and cables to plot
%   cable_format: character string made from elements of plotting columns
%       (e.g. '-k')

% Find vector normal to a ground face
normal_vec = findAvgNormalVector(r, ground_face);
% Want vector to point down, not up (for ground face that's on
% bottom)
if (sign(normal_vec(3)) == 1)
    normal_vec = -normal_vec;
end

% Find rotation matrix to orient ground face downwards
down = [0 0 -1];
normal_unit = normal_vec/norm(normal_vec);
v = cross(down,normal_unit);
s = norm(v);
if s ~= 0
    c = dot(down,normal_unit);
    v_ss = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    R = eye(3) + v_ss + v_ss^2*(1-c)/s^2;
else
    R = eye(3);
end

% Rotate nodes
r_rot = r(:,:,end)*R;

% Find COG of rotated tensegrity
COG = mean(r_rot,1);

% Final tensegrity
figure
plotTensegrity(r_rot, cable_pair, rod_pair, labels_on, cable_format)
hold on
r_rot_ground = r_rot(ground_face,:);
scatter3(r_rot_ground(:,1),r_rot_ground(:,2),r_rot_ground(:,3),'Filled','b')
hold on
scatter3(COG(1),COG(2),COG(3),'Filled',cable_format)
hold on
plot3([COG(1); COG(1)], [COG(2); COG(2)], [COG(3); ...
    min(r(:,3,end))],['--' cable_format])
title(['Cables used: ' num2str(num_cables_used)])