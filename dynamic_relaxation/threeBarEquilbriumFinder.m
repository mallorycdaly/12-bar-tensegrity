%% Equilbrium Finder for Three Bar Tensegrity
% This script finds the equilbrium position of a three bar tensegrity based
% on desired rest lengths. Dynamic relaxation is used to iteratively reach
% the equilbrium configuration from the initial position.
%
% This version of dynamic relaxation models the rods as stiff springs. 
%
% Author: Mallory Daly
% Affiliation: University of California, Berkeley
%              Mechanical Engineering Department
%              NASA Space Technology Research Fellow
% Last Updated: June 21, 2017

clear; close all

%% Design Parameters

% Three-bar tensegrity geometry
edge_length = 1;        % edge length of the top and bottom triangles
height = 1;             % distance between top and bottom
rotation_angle = 45;    % degrees, between top and bottom
rod_radius = 0.1;       % used to check for intersection
[r0, cable_pair, rod_pair, num_nodes, num_cables, num_rods, L0_cable, ...
    L0_rod] = formThreeBar(edge_length, height, rotation_angle);

% Mass and spring constants
m = 10;                             % mass per node
k_rod = 10000;                      % spring constant of the rods
k_spring = 200;                     % spring constant of the springs
L0_spring = zeros(num_cables,1);    % initial length of the springs

% Desired rest lengths: The length of the string in series with the spring
% rest_length = rand(num_cables,1).*L0_cable;
rest_length = ...
   [0.9649;
    0.1576;
    0.9706;
    1.0465;
    0.5307;
    0.8749;
    0.1419;
    0.4218;
    0.9157];

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
intersect_found = zeros(sim_step,1);    % boolean for rod intersection

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
    
    % Check for rod intersection
    intersect_found(i) = checkRodIntersection(r(:,:,i+1), rod_pair, ...
        num_rods, rod_radius);
    
end

% Store updated rod length
L_rod(:,:,end) = norm(r(rod_pair(1,1),:,end) - ...
    r(rod_pair(1,2),:,end));

% Throw warning if a rod intersection was found.
if any(intersect_found) == 1
    warning('Rod intersection was found during the simulation.')
end

%% Output results
% F_total_end = F_total(:,:,end)
% L_rod_end = L_rod(:,:,end)
% r_end = r(:,:,end)
fprintf('\nForce matrix at end of simulation:\n')
disp(F_total(:,:,end))
fprintf('\nNodal positions at end of simulation:\n')
disp(r(:,:,end))
fprintf('\nRod length change at end of simulation:\n')
disp(L_rod(:,:,end)-L_rod(:,:,1))

%% Plot results

% Tensegrity
fig = figure;
fig.OuterPosition = [10 50 750 450];
plotTensegrity(r0, cable_pair, rod_pair, num_nodes, num_cables, ...
    num_rods, labels_on, style_initial)
% addCoordinateSystemToPlot(r, rod_pair, num_rods)
hold on
plotTensegrity(r(:,:,end), cable_pair, rod_pair, num_nodes, num_cables, ...
    num_rods, labels_on, style_equilbrium)
% addForceToPlot(r(:,:,end),F_rod(:,:,end),'r')
% addForceToPlot(r(:,:,end),F_cable(:,:,end),'b')
addForceToPlot(r(:,:,end),F_total(:,:,end),'g')

% Kinetic energy
fig = figure;
fig.OuterPosition = [10 450 750 400];
plot(0:sim_step-1,KE,'LineWidth',1.5);
xlabel('Time step')
ylabel('Kinetic energy')
grid on