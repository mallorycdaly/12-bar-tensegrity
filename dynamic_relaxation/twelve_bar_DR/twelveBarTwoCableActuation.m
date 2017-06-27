%% Two-Cable Actuation Using Base Cable and Secondary Cable
% This script finds the equilbrium position of a 12-bar tensegrity based
% on desired rest lengths, with two cables actuated. One cable is named the
% base cable and is kept constant. A set of secondary cables is defined.
% Each of these cables is paired with the base cable and the equilbrium
% configuration is evaluated for tipping condition, which occurs when the
% COG escapes the base triangle. Dynamic relaxation is used to iteratively
% reach the equilbrium configuration. This version of dynamic relaxation
% models the rods as stiff springs. 
%
% Author: Mallory C. Daly
% Affiliation: University of California, Berkeley
%              Mechanical Engineering Department
%              NASA Space Technology Research Fellow
% Last Updated: June 26, 2017

clear; close all

%% Design Parameters

% Twelve-bar tensegrity cube geometry
rod_radius = 0.01;
scaling_factor = 0.1;

% Actuated cables
actuated_cable_pair = [ 4  9;  % base deformation cable
                        3 21;  % secondary cable 0
                       18 22;  % secondary cable 1
                       16 20;  % secondary cable 2
                        7 23;  % secondary cable 3
                       13  8;  % secondary cable 4
                        6 14;  % secondary cable 5
                        2 10;  % secondary cable 6
                       15 19;  % secondary cable 7
                        0  5;  % secondary cable 8
                        1 12;  % secondary cable 9
                       11 17] + 1;  % secondary cable 10

% Mass and spring constants
m = 10;             % mass per node
k_rod = 1000;       % spring constant of the rods
k_cable = 1000;    % spring constant of the elastic lattice
L0_spring = 0;      % initial length of the springs

% Percent of cable length for calculation of actuated cables' rest lengths
actuation_percent_base = 0.25;
actuation_percent_secondary = 0.25;

% Simulation variables
sim_steps = 1e3;    % length of simulation
del_t = 1e-2;       % time

% Results diplay format
output_results = 0;  % include detailed output of results

% Plotting format for configurations that meet the tipping condition
plot_KE = 0;           % include plots of kinetic energy
plot_initial = 0;      % include plots of initial configurations with final
style_initial = 'y';   % formats plot style of initial tensegrity
style_final = 'r';     % formats plot style of final tensegrity
style_actuated = 'b';  % formats plot style of actuated cables
labels_on = 1;         % include labels of node, cable, and rod numbers

%% Dynamic relaxation (DR)
num_actuated_cables = size(actuated_cable_pair,1);
for i = 1:num_actuated_cables-1

    % Form 12-bar for current actuated cables
    curr_actuated_cable_pair = [actuated_cable_pair(1,:);
                                actuated_cable_pair(i+1,:)];
    [r0, cable_pair, rod_pair, L0_cable, L0_rod, ground_face] = ...
        formTwelveBarCube(scaling_factor, curr_actuated_cable_pair);
    
    % Grab number of unactuated cables
    num_unactuated_cables = size(cable_pair,1) - ...
        size(curr_actuated_cable_pair,1);
    
    % Change length of current actuated cables
    rest_lengths = 0.99*L0_cable;
    rest_lengths(end-1) = actuation_percent_base*L0_cable(end-1);   % base cable
    rest_lengths(end) = actuation_percent_secondary*L0_cable(end);  % secondary cable
    
    % Run DR
    [r, v, KE, F_cable, F_rod, F_total, L_rod] = dynamicRelaxation(r0, ...
        cable_pair, rod_pair, rod_radius, m, k_cable, L0_spring, ...
        k_rod, L0_rod, rest_lengths, sim_steps, del_t);
    
    % Check for rod intersection and throw warning if found
    % [intersect_found, P_intersect, P_distance] = ...
    intersect_found = checkRodIntersection(r(:,:,end), rod_pair, ...
        rod_radius);
    if intersect_found == 1
        fprintf('\n')
        warning(['Rod intersection found in final configuration for ' ...
            'secondary cable ' num2str(i-1)])
    end

    % Output dynamic relaxation results
    if output_results == 1
        fprintf('\nForce matrix at end of simulation:\n')
        disp(F_total(:,:,end))
        fprintf('\nNodal positions at end of simulation:\n')
        disp(r(:,:,end))
        fprintf('\nRod length percent change:\n')
        disp((L_rod(:,:,end)-L_rod(:,:,1))/L0_rod*100)
    end
        
    % Step condition: COG escapes supporting triangle
    num_polygons = size(ground_face,1);
    escaped_poly = zeros(num_polygons,1);
    distance = -inf*ones(num_polygons,1);
    edge_closest = zeros(num_polygons,1);
    for j = 1:num_polygons

        % Find vector normal to a ground face
        normal_vec = findAvgNormalVector(r(:,:,end), ground_face(j,:));

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

        % % Find projected COG
        % COG_proj = [COG(1:2) r_rot(ground_face(i,1),3,end)];

        % Find if projection of COG onto the ground plane is inside or 
        % outside the supporting polygon
        ground_poly = findGroundPoly(r_rot, ground_face(j,:));
        escaped_poly(j) = ~inpolygon(COG(1), COG(2), ground_poly(:,1), ...
            ground_poly(:,2));

        % Find shortest distance to an edge, if it is outside the polygon
        wrapN = @(x,N)(1+mod(x,N));
        if escaped_poly(j) == 1
            [distance(j),~,~,~,edge_closest(j)] = p_poly_dist(COG(1), ...
                COG(2), ground_poly(:,1), ground_poly(:,2), true);
        end

        % Plot tensegrity, if COG escaped
        if escaped_poly(j) == 1
            figure
            plotTensegrity(r(:,:,end), ...
                cable_pair(1:num_unactuated_cables,:), rod_pair, ...
                labels_on, style_final)
            hold on
            plotTensegrity(r(:,:,end), curr_actuated_cable_pair, ...
                rod_pair, 0, style_actuated)
            hold on
            COG = mean(r(:,:,end),1);
            scatter3(COG(1),COG(2),COG(3),'Filled','r')
            title(['Secondary cable: ' num2str(i-1)])
            
            % Initial tensegrity
            if plot_initial == 1
                hold on
                plotTensegrity(r0, cable_pair, rod_pair, 0, style_initial)
            end
            
            % Kinetic energy
            if plot_KE == 1
                figure
                plot(0:sim_steps-1,KE/max(KE),'LineWidth',1.5);
                xlabel('Time step')
                ylabel('Normalized kinetic energy')
                grid on
            end
            
        end
    
    end
    
    % Output roll condition results
    if all(escaped_poly == 0)
        fprintf(['\nThe step condition was NOT met for secondary ' ... 
            'cable %i.\n'], i-1)
    else
        for j = find(escaped_poly == 1)'
            fprintf(['\nThe step condition was met for secondary ' ...
                'cable %i!\nThe distance from the edge was %.5f.\n'], ...
                i-1, distance(j))
        end
    end

end