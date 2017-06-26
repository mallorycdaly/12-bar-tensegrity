%% Equilbrium Finder for 12-Bar Tensegrity
% This script finds the equilbrium position of a 12-bar tensegrity based
% on desired rest lengths. Dynamic relaxation is used to iteratively reach
% the equilbrium configuration from the initial position.
%
% This version of dynamic relaxation models the rods as stiff springs. 
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
                        3 21;
                       18 22;
                       16 20;
                        7 23;
                       13  8;
                        6 14;
                        2 10;
                       15 19;
                        0  5;
                        1 12;
                       11 17] + 1; 

% Mass and spring constants
m = 10;             % mass per node
k_rod = 1000;       % spring constant of the rods
k_lattice = 100;    % spring constant of the elastic lattice
k_cable = 1000;     % spring constant of the actuated cable
L0_spring = 0;      % initial length of the springs

% Percent of cable length for calculation of actuated cables' rest lengths
actuation_percent_base = 0.2;
actuation_percent_secondary = 0.9; %0.9;

% Simulation variables
sim_steps = 1e3;    % length of simulation
del_t = 1e-2;       % time

% Results diplay format
output_results = 0;         % include detailed output of results

% Plotting format
plot_KE = 0;                % include plot of kinetic energy
plot_final = 0;             % plots final tensegrity configuration
plot_initial = 0;           % plots initial configuration with final
plot_rot = 1;               % plots rotated tensegrity for projection of COG
style_initial = 'b';        % formats plot style of initial tensegrity
style_final = 'r';          % formats plot style of final tensegrity
style_actuated = 'b';       % formats plot style of actuated cables
labels_on = 1;              % adds labels of node, cable, and rod numbers

%% Dynamic relaxation (DR)
num_actuated_cables = size(actuated_cable_pair,1);
% for i = 1
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
    rest_lengths(end-1) = actuation_percent_base*L0_cable(end-1);  % base cable
    rest_lengths(end) = actuation_percent_secondary*L0_cable(end);  % secondary cable
    
    % Run DR
    [r, v, KE, F_cable, F_rod, F_total, L_rod] = dynamicRelaxation(r0, ...
        cable_pair, rod_pair, rod_radius, m, k_lattice, L0_spring, k_rod, ...
        L0_rod, rest_lengths, sim_steps, del_t);
    
    % Check for rod intersection and throw warning if found
    % [intersect_found, P_intersect, P_distance] = ...
    intersect_found = checkRodIntersection(r(:,:,end), rod_pair, rod_radius);
    if intersect_found == 1
        fprintf('\n')
        warning(['Rod intersection found in final configuration for secondary cable ' num2str(i)])
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

    % Plot dynamic relaxation results

        % Kinetic energy
        if plot_KE == 1
            figure
            % figure('OuterPosition', [10 500 750 350])
            plot(0:sim_steps-1,KE,'LineWidth',1.5);
            xlabel('Time step')
            ylabel('Kinetic energy')
            grid on
        end

        % Final tensegrity
        if plot_final == 1
            figure        
            plotTensegrity(r(:,:,end), cable_pair(1:num_unactuated_cables,:), ...
                rod_pair, labels_on, style_final)
            hold on
            plotTensegrity(r(:,:,end), curr_actuated_cable_pair, rod_pair, ...
                0, style_actuated)
            addForceToPlot(r(:,:,end),F_total(:,:,end),'g')
            title(['Secondary cable actuated: ' num2str(i)])
        end
        
        % Initial tensegrity
        if plot_initial == 1
            if plot_final == 1
                hold on
            else
                figure
            end
            plotTensegrity(r0, cable_pair, rod_pair, 0, style_initial)
            % addCoordinateSystemToPlot(r, rod_pair, num_rods)  % plot coordinate system
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

        % Find if projection of COG onto the ground plane is inside or outside 
        % the supporting polygon
        ground_poly = findGroundPoly(r_rot, ground_face(j,:));
        escaped_poly(j) = ~inpolygon(COG(1), COG(2), ground_poly(:,1), ...
            ground_poly(:,2));

        % Find shortest distance to an edge, if it is outside the polygon
        wrapN = @(x,N)(1+mod(x,N));
        if escaped_poly(j) == 1
            [distance(j),~,~,~,edge_closest(j)] = p_poly_dist(COG(1), COG(2), ...
                ground_poly(:,1), ground_poly(:,2), true);
        end

        % Plot tensegrity, if COG escaped
        if (escaped_poly(j) == 1 && plot_rot == 1)
            figure
            plotTensegrity(r(:,:,end), cable_pair(1:num_unactuated_cables,:), ...
                rod_pair, labels_on, style_final)
            hold on
            plotTensegrity(r(:,:,end), curr_actuated_cable_pair, rod_pair, ...
                0, style_actuated)
            hold on
            COG = mean(r(:,:,end),1);
            scatter3(COG(1),COG(2),COG(3),'Filled','r')
            hold on
            title(['Secondary cable: ' num2str(i) '. Ground face ' num2str(j) '.'])
        end
    
    end
    
    % Output roll condition results
    if all(escaped_poly == 0)
        fprintf('\nThe step condition was NOT met for secondary cable %i.\n', i)
    else
        for j = find(escaped_poly == 1)'
            fprintf(['\nThe step condition was met for secondary cable %i. ' ...
                '\nThe distance from the edge was %.5f.\n'], ...
                i, distance(j))
        end
    end

end