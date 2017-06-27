%% Equilbrium Finder for 12-Bar Tensegrity
% This script finds the equilbrium position of a 12-bar tensegrity based
% on desired rest lengths. Dynamic relaxation is used to iteratively reach
% the equilbrium configuration from the initial position. This version of
% dynamic relaxation models the rods as stiff springs.
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

% Mass and spring constants
m = 10;             % mass per node
k_rod = 1000;       % spring constant of the rods
k_lattice = 100;    % spring constant of the elastic lattice
k_cable = 1000;     % spring constant of the actuated cable
L0_spring = 0;      % initial length of the springs
[r0, cable_pair, rod_pair, L0_cable, L0_rod, ground_face] = ...
    formTwelveBarCube(scaling_factor, []);

% Rest lengths
rest_lengths = 0.8*L0_cable;

% Simulation variables
sim_steps = 1e3;    % length of simulation
del_t = 1e-2;       % time

% Plotting format
style_initial = 'b';        % formats plot style of initial tensegrity
style_final = 'r';          % formats plot style of final tensegrity
style_actuated = 'k';       % formats plot style of actuated cables
labels_on = 1;              % adds labels of node, cable, and rod numbers
plot_initial = 0;           % plots initial configuration with final

%% Dynamic relaxation (DR)

% Run DR
[r, v, KE, F_cable, F_rod, F_total, L_rod] = dynamicRelaxation(r0, ...
    cable_pair, rod_pair, rod_radius, m, k_lattice, L0_spring, k_rod, ...
    L0_rod, rest_lengths, sim_steps, del_t);

% Output results
fprintf('\nForce matrix at end of simulation:\n')
disp(F_total(:,:,end))
fprintf('\nNodal positions at end of simulation:\n')
disp(r(:,:,end))
fprintf('\nRod length percent change:\n')
disp((L_rod(:,:,end)-L_rod(:,:,1))/L0_rod*100)

% Plot results

    % Kinetic energy
    figure
    plot(0:sim_steps-1,KE,'LineWidth',1.5);
    xlabel('Time step')
    ylabel('Kinetic energy')
    grid on

    % Final tensegrity
    figure
    plotTensegrity(r(:,:,end), cable_pair, rod_pair, ...
        labels_on, style_final)
    addForceToPlot(r(:,:,end),F_total(:,:,end),'g')  % plot total forces
    
    % Initial tensegrity
    if plot_initial == 1
        hold on
        plotTensegrity(r0, cable_pair, rod_pair, 0, style_initial)
        % addCoordinateSystemToPlot(r, rod_pair, num_rods)
        title('Initial and Final Configurations')
    else
        title('Final Configuration')
    end

%% Step condition: COG escapes supporting triangle
% To do: make this generalizable to non-triangle faces
num_polygons = size(ground_face,1);
escaped_poly = zeros(num_polygons,1);
distance = -inf*ones(num_polygons,1);
edge_closest = zeros(num_polygons,1);
for i = 1:num_polygons
    
    % Find vector normal to a ground face
    normal_vec = findAvgNormalVector(r(:,:,end), ground_face(i,:));

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
    ground_poly = findGroundPoly(r_rot, ground_face(i,:));
    escaped_poly(i) = ~inpolygon(COG(1), COG(2), ground_poly(:,1), ...
        ground_poly(:,2));

    % Find shortest distance to an edge, if it is outside the polygon
    wrapN = @(x,N)(1+mod(x,N));
    if escaped_poly(i) == 1
        [distance(i),~,~,~,edge_closest(i)] = p_poly_dist(COG(1), COG(2), ...
            ground_poly(:,1), ground_poly(:,2), true);
    end
    
    % Plot COG
    hold on
    scatter3(COG(1),COG(2),COG(3),'Filled','r')
    
    % Plot rotated tensegrity, if COG escaped
    if escaped_poly(i) == 1
        figure
    %     figure('OuterPosition',[750 50 750 450])
    %     plotTensegrity(r0, cable_pair, rod_pair, num_nodes, num_cables, ...
    %         num_rods, labels_on, style_initial)
    %     hold on
        plotTensegrity(r_rot, cable_pair, rod_pair, labels_on, ...
            style_equilbrium)
        hold on
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