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

% clear; close all
clear

%% Design Parameters

% Twelve-bar tensegrity cube geometry
scaling_factor = 0.1;   % scales node positions
rod_radius = 0.05;      % used for rod intersection check, approximately 3/16 in.
% ground_face = [4 6 3 23 9 18 21 15] + 1;    % octagon base of cube
ground_face = [9 20 23] + 1;                % triangle base of cube
% ground_face = [4 5 6 12 17 11] + 1;         % hexagon base of octahedron
% ground_face = [17 12 18 23] + 1;            % square base of octahedron
cross_body_pair = [];   % node index pairs that define cross-body cables
[r0, cable_pair, rod_pair, L0_cable, L0_rod] = formTwelveBarCube(...
    cross_body_pair);
% [r0, cable_pair, rod_pair, L0_cable, L0_rod] = formTwelveBarOctahedron(...
%     cross_body_pair);

% Mass and spring constants
m = 1;              % mass per node (tuned for convergence)
k_rod = 4000;       % spring constant of the rods (tuned for convergence)
k_cable = 500;      % spring constant of the elastic lattice (tuned for convergence)
L0_spring = 0;      % initial length of the springs (careful: no checks if L0_spring is too long)

% Percent of initial cable length for rest length of actuated cables
percent_length = 0.05;

% Number of actuated cables
num_act_cables = 1;

% Simulation variables
sim_steps = 1e3;    % length of simulation
del_t = 1e-2;       % time

% Results diplay format
output_results = 1;    % include detailed output of results

% Plotting format for configurations that meet the tipping condition
plot_initial = 0;      % include plots of initial configurations with final
labels_on = 1;         % include labels of node, cable, and rod numbers
color_initial = 'g';   % formats plot color of initial tensegrity
color_final = 'r';     % formats plot color of final tensegrity
color_actuated = 'b';  % formats plot color of actuated cables
plot_KE = 1;           % include plots of kinetic energy

%% Dynamic relaxation (DR)

% Find cable combinations and grab numbers
num_cables = size(cable_pair,1);
cable_combos = combnk(1:num_cables,num_act_cables);
num_combos = size(cable_combos,1);
num_rods = size(rod_pair,1);
num_nodes = size(r0,1);

% Initialize results storage
results.cable_combos = cable_combos-1;
results.edge_closest = zeros(num_combos,1);
results.edge_distance = zeros(num_combos,1);
results.ground_face = ground_face-1;
results.escaped = zeros(num_combos,1);
results.intersect = zeros(num_combos,1);
results.r = zeros([size(r0) num_combos]);           % 3D array
results.r_rot = zeros([size(r0) num_combos]);       % 3D array
results.F_cable = zeros([size(r0) num_combos]);     % 3D array
results.F_rod = zeros([size(r0) num_combos]);       % 3D array
results.KE = zeros(sim_steps,1,num_combos);         % 3D array
results.L_rod = zeros(num_rods,1,num_combos);       % 3D array

% Begin dynamic relaxation
for i = 1:size(cable_combos,1)
    
    % Set rest lengths
    rest_lengths = 0.95*L0_cable;
    rest_lengths(cable_combos(i,:)) = percent_length * ...
        L0_cable(cable_combos(i,:));
    
    % Run dynamic relaxation
    [r, ~, KE, F_cable, F_rod, ~, L_cable, L_rod, intersect_found] = ...
        dynamicRelaxation(r0, cable_pair, rod_pair, m, k_cable, ...
        L0_spring, k_rod, L0_rod, rod_radius, rest_lengths, sim_steps, ...
        del_t);

    % Check tipping condition: COG escapes supporting polygon
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
        results.edge_closest(i) = edge_closest;
        results.edge_distance(i) = distance;
        results.escaped(i) = escaped_poly;
        results.intersect(i) = intersect_found;
        results.r(:,:,i) = r(:,:,end);
        results.r_rot(:,:,i) = r_rot;
        results.F_cable(:,:,i) = F_cable(:,:,end);
        results.F_rod(:,:,i) = F_rod(:,:,end);
        results.KE(:,:,i) = KE;
        results.L_cable(:,:,i) = L_cable(:,:,end);
        results.L_rod(:,:,i) = L_rod(:,:,end);
        
end

%% Analyze results
% Goal to find cable combinations to move over each edge of the supporting
% polygon. This section outputs those cable combinations that resulted in
% the most displacement of the COG from each edge of the supporting
% polygon.
close all

% Find indices corresponding to max distance from each edge
max_idx = findMaxIdxEachEdge(results);

% Plot results
for i = 1:length(max_idx)
    if max_idx(i) ~= 0
        
        % Final tensegrity (convert to cm)
        figure
        plotTensegrity(results.r_rot(:,:,max_idx(i))*10, cable_pair, ...
            rod_pair, labels_on, color_final)
%         addForceToPlot(results.r_rot(:,:,max_idx(i)), ...
%             results.F_rod(:,:,max_idx(i)) + ...
%             results.F_cable(:,:,max_idx(i)),'g')  % plot total forces
        % Add COG
        COG = mean(results.r_rot(:,:,max_idx(i)),1)*10;
        hold on
        scatter3(COG(1),COG(2),COG(3),'Filled',color_final)
        % Add COG's projection to ground
        hold on
        plot3([COG(1); COG(1)], [COG(2); COG(2)], [COG(3); ...
            min(results.r_rot(:,3,max_idx(i)))*10],['--' color_final])
        % Highlight actuated cable(s)
        curr_cable_combo = results.cable_combos(max_idx(i),:) + 1;
        for j = 1:length(curr_cable_combo)
            hold on
            act_cable_pair = cable_pair(curr_cable_combo(j),:);
            plotTensegrity(results.r_rot(:,:,max_idx(i))*10, ...
                act_cable_pair, rod_pair, 0, color_actuated)
        end
        xlabel('x (cm)')
        ylabel('y (cm)')
        zlabel('z (cm)')
%         set(gca,'fontsize',12)
        grid on
        
        % Kinetic energy
        if plot_KE == 1
            figure
            curr_KE = results.KE(:,:,max_idx(i));
%             plot(0:sim_steps-1,curr_KE/max(curr_KE),'LineWidth',1.5);
            plot(0:300-1,curr_KE(1:300)/max(curr_KE),'LineWidth',1.25);
            xlabel('Time step')
            ylabel('Normalized kinetic energy')
%             set(gca,'fontsize',12)
            grid on
        end
        
        % Output results
        if output_results == 1
            
            fprintf(['\n-----------------------------------------------'...
                '------------------'])
            
            fprintf('\nActuated cable(s):\n')
            disp(results.cable_combos(max_idx(i),:))
            
            fprintf('\nDistance COG escaped (cm):\n')
            disp(results.edge_distance(max_idx(i))*10)

            fprintf('\nMaximum final total force (N):\n')
            disp(max(max((results.F_rod(:,:,max_idx(i)) + ...
                results.F_cable(:,:,max_idx(i)))))/10)
            
            fprintf('\nFinal KE (J):\n')
            disp(curr_KE(end)/100)
            
            fprintf('\nMaximum rod length percent change:\n')
            disp(min((results.L_rod(:,:,max_idx(i))-L0_rod)/L0_rod*100))
            
            if results.intersect(max_idx(i)) == 1
                fprintf('\nWarning: intersection found on this edge.\n')
            end
            
        end
        
    end
end
