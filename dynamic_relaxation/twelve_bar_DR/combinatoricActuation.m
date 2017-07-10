%% Combinatoric Approach to Finding Actuation Scheme
% This script finds the equilbrium positions of a 12-bar tensegrity for a
% set of actuated cables. The set of cables is varied systematically using
% combinatorics. The user inputs the rest length of the actuated cables
% and the number of cables to be actuated. The program finds the equilbrium
% positions that result from all combinations of cables available. The
% tipping condition is evaluated for every combination and the results are
% stored.
%
% Dynamic relaxation is used to iteratively reach equilbrium. This version
% of dynamic relaxation models the rods as stiff springs. 
%
% Note: Node, cable, and rod numbers use zero indexing for easy transition
% to NTRT. Thus "+1" is added to convert to MATLAB indexing.
%
% Author: Mallory C. Daly
% Affiliation: University of California, Berkeley
%              Mechanical Engineering Department
%              NASA Space Technology Research Fellow
% Created: June 2017
% Last Updated: July 2017

clear; close all

%% Design Parameters

% Twelve-bar tensegrity cube geometry
scaling_factor = 0.1;   % scales node positions
rod_radius = 0.01;      % used for rod intersection check
ground_face = [4 6 3 23 9 18 21 15] + 1;    % octagon base of cube
% ground_face = [9 20 23] + 1;                % triangle base of cube
% ground_face = [4 5 6 12 17 11] + 1;         % hexagon base of octahedron
cross_body_pair = [];   % node index pairs that define cross-body cables
[r0, cable_pair, rod_pair, L0_cable, L0_rod] = formTwelveBarCube(...
    cross_body_pair);
% [r0, cable_pair, rod_pair, L0_cable, L0_rod] = formTwelveBarOctahedron(...
%     cross_body_pair);

% Mass and spring constants
m = 10;             % mass per node
k_rod = 1000;       % spring constant of the rods
k_cable = 200;      % spring constant of the elastic lattice
L0_spring = 0;      % initial length of the springs

% Percent of initial cable length for rest length of actuated cables
percent_length = 0.2;

% Combinatoric variables
min_choose_k = 1;   % min number of cables pulled at once
max_choose_k = 3;   % max number of cables pulled at once

% Simulation variables
sim_steps = 1e3;    % length of simulation
del_t = 1e-2;       % time

% Results diplay format
output_results = 0;  % include detailed output of results

% Plotting format for configurations that meet the tipping condition
plot_KE = 0;           % include plots of kinetic energy
plot_initial = 0;      % include plots of initial configurations with final
color_initial = 'g';   % formats plot color of initial tensegrity
color_final = 'r';     % formats plot color of final tensegrity
color_actuated = 'b';  % formats plot color of actuated cables
labels_on = 1;         % include labels of node, cable, and rod numbers

%% Dynamic relaxation (DR)

% Loop through cable combinations
num_cables = size(cable_pair,1);
results = {};
for i = min_choose_k:max_choose_k
    cable_combos = combnk(1:num_cables,i);
    results.escaped{i} = zeros(size(cable_combos,1),3);
    results.r{i} = zeros(size(r0,1),3,size(cable_combos,1));

for j = 1:size(cable_combos,1)
    
    % Set rest lengths
    rest_lengths = 0.95*L0_cable;
    rest_lengths(cable_combos(j,:)) = percent_length * ...
        L0_cable(cable_combos(j,:));
    
    % Run dynamic relaxation
    [r, ~, KE, F_cable, F_rod, F_total, L_rod, intersect_found] = ...
        dynamicRelaxation(r0, cable_pair, rod_pair, m, k_cable, ...
        L0_spring, k_rod, L0_rod, rod_radius, rest_lengths, sim_steps, ...
        del_t);

    % Step condition: COG escapes supporting polygon
    distance = -inf;
    edge_closest = 0;

        % Find vector normal to a ground face
        normal_vec = findAvgNormalVector(r(:,:,end), ground_face);
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

        % Find if projection of COG onto the ground plane is inside or 
        % outside the supporting polygon
        ground_poly = r_rot(ground_face,:);
        escaped_poly = ~inpolygon(COG(1), COG(2), ground_poly(:,1), ...
            ground_poly(:,2));

        % Find shortest distance to an edge, if it is outside the polygon
        wrapN = @(x,N)(1+mod(x,N));
        if escaped_poly == 1
            [distance,~,~,~,edge_closest] = p_poly_dist(COG(1), ...
                COG(2), ground_poly(:,1), ground_poly(:,2), true);
        end
        
        % Store results
        results.escaped{i}(j,:) = [escaped_poly distance edge_closest];
        results.r{i}(:,:,j) = r(:,:,end);
        results.r_rot{i}(:,:,j) = r_rot;
        results.F_cable{i}(:,:,j) = F_cable(:,:,end);
        results.F_rod{i}(:,:,j) = F_rod(:,:,end);
        results.KE{i}(:,:,j) = KE;
        results.intersect{i}(j,:) = intersect_found;
        
end

        % Store cable combos
        results.cable_combos{i} = cable_combos-1;

end

%% Analyze results
% TO DO: Make this only plot max distance result
% Plot results that escaped base polygon
% color = 'r';
% for i = min_choose_k:max_choose_k
%     postAnalyzer(results, i, ground_face, cable_pair, rod_pair, ...
%         labels_on, color)
% end