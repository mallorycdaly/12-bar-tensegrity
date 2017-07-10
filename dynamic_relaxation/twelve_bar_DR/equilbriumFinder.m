%% 12-Bar Equilbrium Finder
% This script finds the equilbrium position of a 12-bar tensegrity based
% on desired rest lengths. Dynamic relaxation is used to iteratively reach
% the equilbrium configuration from the initial position. This version of
% dynamic relaxation models the rods as stiff springs.
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
rod_radius = 0.01;      % used for rod intersection check
scaling_factor = 0.1;   % scales node positions

% Mass and spring constants
m = 10;             % mass per node
k_rod = 1000;       % spring constant of the rods
k_lattice = 200;    % spring constant of the elastic lattice
L0_spring = 0;      % initial length of the springs

% Stored info
% 
% 12-bar cube:
% ground_face = [ 1 16  5 18  9 20 13 22;
%                 0  2 15 21  5 16 19 11;
%                 0  2  4  6  8 10 12 14;
%                10  8  3 23 20 13 17  7;
%                14 12  7 17 22  1 19 11;
%                 4  6  3 23  9 18 21 15] + 1;
% cross_body_pair = [ 4  9;         % bottom /
%                     1 12;         % top /
%                    18 22;         % front /
%                    15 19;         % left /
%                     6 14;         % back /
%                     7 23;         % right /
%                     3 21;         % bottom \
%                    11 17;         % top \
%                    16 20;         % front \
%                     0  5;         % left \
%                     2 10;         % back \
%                     8 13] + 1;  	% right\

% Form 12-bar
cross_body_pair = [];
[r0, cable_pair, rod_pair, L0_cable, L0_rod] = formTwelveBarOctahedron(...
    cross_body_pair);

% Rest lengths
% percent_length = rand(size(L0_cable))
% rest_lengths = percent_length.*L0_cable;
rest_lengths = 0.95*L0_cable;
% rest_lengths([1 12 35 36 28 16 20 24 32 7 11 33]) = 0.3*L0_cable([1 12 35
% 36 28 16 20 24 32 7 11 33]);  % make triangles larger
% rest_lengths(end) = 0.2*L0_cable(end);  % actuate last segment

% Simulation variables
sim_steps = 1e3;    % length of simulation
del_t = 1e-2;       % time

% Plotting format
color_initial = 'b';        % formats color of initial tensegrity plot
color_final = 'r';          % formats color of final tensegrity plot
color_actuated = 'k';       % formats color of actuated cables
labels_on = 0;              % add labels of node, cable, and rod numbers
plot_initial = 0;           % plot initial configuration with final
plot_KE = 1;                % plot kinetic energy

%% Dynamic relaxation (DR)

% Run DR
[r, v, KE, F_cable, F_rod, F_total, L_rod, intersect_found] = ...
    dynamicRelaxation(r0, cable_pair, rod_pair, m, k_lattice, ...
    L0_spring, k_rod, L0_rod, rod_radius, rest_lengths, sim_steps, del_t);

% Throw warning if rod intersection was found
if intersect_found == 1
    fprintf('\n')
    warning(['Rod intersection found in final configuration for ' ...
        'secondary cable ' num2str(i-1)])
end

% Output results
fprintf('\nForce matrix at end of simulation:\n')
disp(F_total(:,:,end))
% fprintf('\nNodal positions at end of simulation:\n')
% disp(r(:,:,end))
fprintf('\nRod length percent change:\n')
disp((L_rod(:,:,end)-L_rod(:,:,1))/L0_rod*100)

% Plot results

    % Kinetic energy
    if plot_KE == 1
        figure
        plot(0:sim_steps-1,KE/max(KE),'LineWidth',1.5);
        xlabel('Time step')
        ylabel('Normalized kinetic energy')
        grid on
    end

    % Final tensegrity
    figure
    plotTensegrity(r(:,:,end), cable_pair, rod_pair, labels_on, ...
        color_final)
    addForceToPlot(r(:,:,end),F_total(:,:,end),'g')  % plot total forces
    COG = mean(r(:,:,end),1);
    hold on
    scatter3(COG(1),COG(2),COG(3),'Filled',color_final)
    hold on
    plot3([COG(1); COG(1)], [COG(2); COG(2)], [COG(3); ...
        min(r(:,3,end))],['--' color_final])
    
    % Initial tensegrity
    if plot_initial == 1
        hold on
        plotTensegrity(r0, cable_pair, rod_pair, labels_on, color_initial)
        % addCoordinateSystemToPlot(r, rod_pair, num_rods)
        title('Initial and Final Configurations')
    else
        title('Final Configuration')
    end

%% Step condition: COG escapes supporting polygon
num_polygons = size(ground_face,1);
escaped_poly = zeros(num_polygons,1);
distance = -inf*ones(num_polygons,1);
edge_closest = zeros(num_polygons,1);
for i = 1:num_polygons
    
    % Find vector normal to a ground face
    normal_vec = findAvgNormalVector(r(:,:,end), ground_face(i,:));
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
    if round(s,12) ~= 0
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
    ground_poly = findGroundPoly(r_rot, ground_face(i,:));
    escaped_poly(i) = ~inpolygon(COG(1), COG(2), ground_poly(:,1), ...
        ground_poly(:,2));

    % Find shortest distance to an edge, if it is outside the polygon
    wrapN = @(x,N)(1+mod(x,N));
    if escaped_poly(i) == 1
        [distance(i),~,~,~,edge_closest(i)] = p_poly_dist(COG(1), COG(2), ...
            ground_poly(:,1), ground_poly(:,2), true);
    end
    
    % Plot rotated tensegrity, if COG escaped
    if escaped_poly(i) == 1
        figure
        plotTensegrity(r_rot, cable_pair, rod_pair, labels_on, ...
            color_final)
        hold on
        scatter3(COG(1),COG(2),COG(3),'Filled',color_final)
        title(['Ground face ' num2str(i)])
    end
    
end

%% Output roll condition results
if all(escaped_poly == 0)
    fprintf('\nThe step condition was NOT met.\n')
else
    for i = find(escaped_poly == 1)'
        fprintf(['\nThe step condition was met for a projection onto ' ...
            'ground face %i!\nThe distance from the edge was %.5f.\n'], ...
            i, distance(i))
    end
end