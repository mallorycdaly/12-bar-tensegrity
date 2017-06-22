%% Equilbrium Finder for 12-Bar Tensegrity
% This script finds the equilbrium position of a 12-bar tensegrity based
% on desired rest lengths. Dynamic relaxation is used to iteratively reach
% the equilbrium configuration from the initial position.
%
% This version of dynamic relaxation models the rods as stiff springs. 
%
% Author: Mallory Daly
% Affiliation: University of California, Berkeley
%              Mechanical Engineering Department
%              NASA Space Technology Research Fellow
% Last Updated: June 22, 2017

clear; close all

%% Design Parameters

% Twelve-bar tensegrity cube geometry
rod_radius = 0.01;
scaling_factor = 0.1;
[r0, cable_pair, rod_pair, num_nodes, num_cables, num_rods, L0_cable, ...
    L0_rod, ground_face] = formTwelveBarCube(scaling_factor);

% Mass and spring constants
m = 3;                             % mass per node
k_rod = 1000;                      % spring constant of the rods
k_spring = 100;                     % spring constant of the springs
L0_spring = zeros(num_cables,1);    % initial length of the springs


% Desired rest lengths: The length of the string in series with the spring
rest_length = 0.99*L0_cable;
rest_length(end) = 0.05*L0_cable(end);
% rest_length([36 29 34 35]) = 0.1*L0_cable([36 29 34 35]);

% rest_length = rand(num_cables,1).*L0_cable;

% Simulation variables
sim_step = 1e3;     % length of simulation
del_t = 1e-2;       % time

% Plotting format
style_initial = 'b';        % formats plot style of initial tensegrity
style_equilbrium = 'r';     % formats plot style of equilbrium tensegrity
labels_on = 1;              % adds labels of node, cable, and rod numbers

%% Dynamic relaxation

% Initialize variables
% Note: All variables except velocity are incremented by del_t starting
% from t = 0 (index 1 in MATLAB). Velocity is incremented by del_t starting
% from t = -del_t/2 (index 1).
r = zeros(num_nodes,3,sim_step+1);      % nodal positions
r(:,:,1) = r0;                          % initial nodal positions are known
v = zeros(num_nodes,3,sim_step+1);    	% nodal velocities
KE = zeros(sim_step,1);                 % kinetic energy
F_cable = zeros(num_nodes,3,sim_step);  % cable force at nodes
F_rod = zeros(num_nodes,3,sim_step);    % rod force at nodes
F_total = zeros(num_nodes,3,sim_step);  % total force at nodes
L_rod = zeros(num_rods,1,sim_step+1);   % rod length

% Run dynamic relaxation
for i = 1:sim_step

    % Find forces acting on each node
    F_cable(:,:,i) = findCableForce(r(:,:,i), cable_pair, num_nodes, ...
        k_spring, L0_spring, rest_length);
    [F_rod(:,:,i),L_rod(:,:,i)] = findRodForce(r(:,:,i), rod_pair, ...
        num_nodes, num_rods, k_rod, L0_rod);
    F_total(:,:,i) = F_cable(:,:,i) + F_rod(:,:,i);
    
    % Find velocity at t + del_t/2, and velocity centered at t using
    % average of t+del_t/2 and t-del_t/2
    v(:,:,i+1) = v(:,:,i) + F_total(:,:,i)/(m/del_t);
    if i == 1
        v_center = zeros(size(v(:,:,i)));
    else
        v_center = (v(:,:,i+1) + v(:,:,i)) / 2;
    end

    % Find kinetic energy (KE)
    for j = 1:num_nodes
        KE(i) = KE(i) + 0.5*m*v_center(j,:)*v_center(j,:)';
    end
    
    % Kinetic damping: Reset velocity and KE to zero if peak was detected
    % and flag for form finding process to restart
    if i > 1
        if KE(i) - KE(i-1) < 0
            % Reset velocity and KE
            v(:,:,i+1) = 0;
            KE(i) = 0;
        end
    end
    
    % Update position
    r(:,:,i+1) = r(:,:,i) + v(:,:,i+1)*del_t;
    
end

% Store updated rod length
L_rod(:,:,end) = norm(r(rod_pair(1,1),:,end) - ...
    r(rod_pair(1,2),:,end));

% Check for rod intersection
[intersect_found, P_intersect, P_distance] = ...
    checkRodIntersection(r(:,:,end), rod_pair, num_rods, rod_radius);

% Throw warning if a rod intersection was found.
if intersect_found == 1
    warning('Configuration state has intersecting rods.')
end

% %% Output dynamic relaxation results
fprintf('\nForce matrix at end of simulation:\n')
disp(F_total(:,:,end))
fprintf('\nNodal positions at end of simulation:\n')
disp(r(:,:,end))
fprintf('\nRod length percent change:\n')
disp((L_rod(:,:,end)-L_rod(:,:,1))/L0_rod*100)

%% Plot dynamic relaxation results

% Kinetic energy
figure
% figure('OuterPosition', [10 500 750 350])
plot(0:sim_step-1,KE,'LineWidth',1.5);
xlabel('Time step')
ylabel('Kinetic energy')
grid on

% Initial and final tensegrity configurations
figure
% figure('OuterPosition', [10 50 750 450])
% plotTensegrity(r0, cable_pair, rod_pair, num_nodes, num_cables, ...
%     num_rods, labels_on, style_initial)
% addCoordinateSystemToPlot(r, rod_pair, num_rods)  % plot coordinate system
% hold on
plotTensegrity(r(:,:,end), cable_pair, rod_pair, num_nodes, num_cables, ...
    num_rods, labels_on, style_equilbrium)
% addForceToPlot(r(:,:,end),F_rod(:,:,end),'r')       % plot rod forces
% addForceToPlot(r(:,:,end),F_cable(:,:,end),'b')     % plot cable forces
addForceToPlot(r(:,:,end),F_total(:,:,end),'g')     % plot total forces
title('Initial (Blue) and Final (Red) Configurations')

%% Step condition: COG escapes supporting triangle
% To do: make this generalizable to non-triangle faces
escaped_poly = zeros(size(ground_face,1),1);
distance = -inf*ones(size(ground_face,1),1);
edge_closest = zeros(size(ground_face,1),1);
for i = 1:size(ground_face,1)
    
    % Find vector normal to a ground face
    nodeA = r(ground_face(i,1),:,end);
    nodeB = r(ground_face(i,2),:,end);
    nodeC = r(ground_face(i,3),:,end);
    normal_vec = cross(nodeA-nodeB,nodeC-nodeB);

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

    % Find if projection of COG onto the ground plane is inside or outside the 
    % supporting polygon
    ground_poly = [r_rot(ground_face(i,1),:,end);
                   r_rot(ground_face(i,2),:,end);
                   r_rot(ground_face(i,3),:,end)];
    escaped_poly(i) = ~inpolygon(COG(1), COG(2), ground_poly(:,1), ...
        ground_poly(:,2));

    % Find shortest distance to an edge, if it is outside the polygon
    edges = size(ground_face,2);
    wrapN = @(x,N)(1+mod(x,N));
    if escaped_poly(i) == 1
        [distance(i),~,~,~,edge_closest(i)] = p_poly_dist(COG(1), COG(2), ...
            ground_poly(:,1), ground_poly(:,2), true);
    end
    
    % Plot rotated tensegrity
    figure
%     figure('OuterPosition',[750 50 750 450])
    plotTensegrity(r0, cable_pair, rod_pair, num_nodes, num_cables, ...
        num_rods, labels_on, style_initial)
    hold on
    plotTensegrity(r_rot, cable_pair, rod_pair, num_nodes, num_cables, ...
        num_rods, labels_on, style_equilbrium)
    hold on
    scatter3(COG(1),COG(2),COG(3),'Filled','r')
    title(['Ground face ' num2str(i)])
    
end

%% Output roll condition results
if all(escaped_poly == 0)
    fprintf('\nThe step condition was NOT met.\n')
else
    for i = find(escaped_poly == 1)'
        fprintf(['\nThe step condition was met for a projection onto ' ...
            'plane %i!\nThe distance from the edge was %1.5f.\n'], i, ...
            distance(i))
    end
end