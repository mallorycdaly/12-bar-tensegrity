%% Test optimization-based equilbrium finder using three strut prism
clear; close all

%% TO DO:
% Add definitions section?
% Add quick access design parameters section

%% Set up three-strut prism tensegrity
% Using prism numbering scheme from Fig. 10 of Review of Form-Finding
% Methods for Tensegrity Structures by A. Tibert.

% Assign design parameters
edge_length = 1;
height = 1;
rotation_angle = 35;  % degrees

% Find initial position of nodes (r0) and indices corresponding to cables
% (cable_pairs) and rods (rod_pairs)
[r0, cable_pair, rod_pair, num_nodes, num_cables, num_rods] = ...
    formThreeBar(edge_length, height, rotation_angle);

% Original lengths, before any deformation
L0_cable = zeros(num_cables,1);
for i = 1: num_cables
    L0_cable(i) = norm(r0(cable_pair(i,1),:) - ...
        r0(cable_pair(i,2),:));
end
L0_rod = norm(r0(rod_pair(1,1),:) - r0(rod_pair(1,2),:));

%% Plot original tensegrity
% figure
% labels_on = 1;
% plotTensegrity(r0, cable_pair, rod_pair, num_nodes, num_cables, ...
%     num_rods, labels_on)
% addCoordinateSystemToPlot()

%% Set up dynamic relaxation

% Set simulation variables
sim_step = 1e3;
del_t = 1e-2;

% Initialize variables
% Note: All variables except velocity are incremented by del_t starting
% from t=0 (index 1 in MATLAB). Velocity is incremented by del_t starting
% from t=del_t/2 (index 1).
r = zeros(num_nodes,3,sim_step+1);      % nodal positions
r(:,:,1) = r0;                          % initial nodal positions are known
v = zeros(num_nodes,3,sim_step);    	% nodal velocities
KE = zeros(sim_step,1);                 % kinetic energy
F_cable = zeros(num_nodes,3,sim_step);  % cable force at nodes
F_rod = zeros(num_nodes,3,sim_step);    % rod force at nodes
F_total = zeros(num_nodes,3,sim_step);  % total force at nodes
L_rod = zeros(num_rods,1,sim_step+1);   % rod length

% Set mass and spring constants
m = 10;                                 % mass per node
k_rod = 1000;                           % spring constant of the rods
k_spring = 10;                          % spring constant of the springs
L0_spring = zeros(num_cables,1);        % initial length of the springs
% c = 1;                                  % damping constant

% Set desired rest lengths
% rest_length = rand(num_cables,1).*L0_cable;
% rest_length = 0.9*L0_cable;
rest_length = 0.75*[1 1 1 0.5 0.5 0.5 1 1 1]'.*L0_cable;

% Set contact force
g = 9.81;
F_contact = [zeros(3,3); [zeros(3,2) 6*m/3*ones(3,1)]];
F_weight = -[zeros(num_nodes,2) m*ones(num_nodes,1)];

%% Dynamic relaxation

% restart = 0;
for i = 1:sim_step

    % Find forces acting on each node
    F_cable(:,:,i) = findCableForce(r(:,:,i), cable_pair, num_nodes, ...
        k_spring, L0_spring, rest_length);
    [F_rod(:,:,i),L_rod(:,:,i)] = findRodForce(r(:,:,i), rod_pair, ...
        num_nodes, num_rods, k_rod, L0_rod);
%     F_total(:,:,i) = F_cable(:,:,i) + F_rod(:,:,i) + F_contact + F_weight;
    F_total(:,:,i) = F_cable(:,:,i) + F_rod(:,:,i);
    
    % Find velocity at t + del_t/2, and velocity centered at t using
    % average of t+del_t/2 and t-del_t/2
    
        % First or restarted velocity update is modified
%         if (i == 1 || restart == 1)
        if i == 1
            v(:,:,i) = F_total(:,:,i)/(m/del_t);            
%             v(:,:,i) = F_total(:,:,i)/(m/del_t+c/2);
            v_center = zeros(size(v(:,:,1)));
            restart = 0;
        else
            v(:,:,i) = v(:,:,i-1) + F_total(:,:,i)/(m/del_t);
%             v(:,:,i) = v(:,:,i-1)*(m/del_t-c/2)/(m/del_t+c/2) + ...
%                 F_total(:,:,i)/(m/del_t+c/2);
            v_center = (v(:,:,i) + v(:,:,i-1)) / 2;
        end
    
    % Find kinetic energy (KE)
    for j = 1:num_nodes
        KE(i) = KE(i) + 0.5*m*v_center(j,:)*v_center(j,:)';
    end
    
    % Kinetic damping: Reset velocity to zero if peak was detected and flag
    % for form finding process to restart
    if i > 1
        if KE(i) - KE(i-1) < 0
            % Reset velocity and restart form finding
            v(:,:,i) = 0;
            KE(i) = 0;
            restart = 1;
%             i
        end
    end
    
    % Update position
    r(:,:,i+1) = r(:,:,i) + v(:,:,i)*del_t;
    
end

% Store updated rod length
L_rod(:,:,end) = norm(r(rod_pair(1,1),:,end) - ...
    r(rod_pair(1,2),:,end));

%% Output results
% F_total(:,:,end)
% L_rod(:,:,end)
% r(:,:,end)

%% Plot results

% Tensegrity
fig = figure;
fig.OuterPosition = [10 50 750 450];
labels_on = 0;
plotTensegrity(r0, cable_pair, rod_pair, num_nodes, num_cables, ...
    num_rods, labels_on, 'b')
% addCoordinateSystemToPlot(r, rod_pair, num_rods)
hold on
plotTensegrity(r(:,:,end), cable_pair, rod_pair, num_nodes, num_cables, ...
    num_rods, labels_on, 'r')
% addForceToPlot(r(:,:,end),F_rod(:,:,end),'r')
% addForceToPlot(r(:,:,end),F_cable(:,:,end),'b')
addForceToPlot(r(:,:,end),F_total(:,:,end),'g')

% Kinetic energy
fig = figure;
fig.OuterPosition = [10 450 750 400];
plot(0:sim_step-1,KE);
xlabel('Time step')
ylabel('Kinetic energy')
grid on

% % Velocity
% norm_v = zeros(sim_step,1);
% for i = 1:sim_step
%     norm_v(i) = norm(v(:,:,i));
% end
% figure
% plot(1:sim_step,norm_v)