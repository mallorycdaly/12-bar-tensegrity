function plotEscaped(results, ground_face, cable_pair, rod_pair, labels_on, color)
% Plots the results that escaped the base polygon.
%
% Inputs:
%   results: structure with two fields, "escaped" and "r"
%   num_comb: number of cables used in combination
%   ground_face: node indices corresponding to ground face
%   cable_pair: node indices corresponding to each cable
%   rod_pair: node indices corresponding to each rod
%   labels_on: boolean to add node/rod/cable labels to plot
%   color: color of cables on plot

close all

% Find results that escaped the base polygon
row = find(results.escaped{num_comb}(:,1));

for i = 1:length(row)
    
    % Grab nodal positions from result
    r = results.r{num_comb}(:,:,row(i));

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
    plotTensegrity(r_rot, cable_pair, rod_pair, labels_on, color)
    % addForceToPlot(r(:,:,end),F_total(:,:,end),'g')  % plot total forces
    hold on
    r_rot_ground = r_rot(ground_face,:);
    scatter3(r_rot_ground(:,1),r_rot_ground(:,2),r_rot_ground(:,3),'Filled','b')
    hold on
    scatter3(COG(1),COG(2),COG(3),'Filled',color)
    hold on
    plot3([COG(1); COG(1)], [COG(2); COG(2)], [COG(3); ...
        min(r_rot(:,3))],['--' color])
    title(['Cables used: ' num2str(num_comb) ', Combination #: ' num2str(row(i))])
    
end