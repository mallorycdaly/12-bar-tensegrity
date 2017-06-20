%% Test optimization-based equilbrium finder using three strut prism
clear; close all

%% Set up three-strut prism tensegrity
% Using prism numbering scheme from Fig. 10 of Review of Form-Finding
% Methods for Tensegrity Structures by A. Tibert.

% Assign design parameters
edge_length = 1;
height = 1;
rotation_angle = -25;  % degrees

% Find initial position of nodes (r0) and indices corresponding to cables
% (cable_pairs) and rods (rod_pairs)
[r0, cable_pair, rod_pair, num_nodes, num_cables, num_rods] = ...
    formThreeBar(edge_length, height, rotation_angle);

% Original lengths, before any deformation
original_cable_length = zeros(num_cables,1);
for i = 1: num_cables
    original_cable_length(i) = norm(r0(cable_pair(i,1),:) - ...
        r0(cable_pair(i,2),:));
end
original_rod_length = norm(r0(rod_pair(1,1),:) - r0(rod_pair(1,2),:));

%% Plot original tensegrity
% figure
% labels_on = 1;
% plotTensegrity(r0, cable_pair, rod_pair, num_nodes, num_cables, ...
%     num_rods, labels_on)
% addCoordinateSystemToPlot()

%% Set up dynamic relaxation

% Set simulation variables
sim_step = 1e5;
del_t = 1e-2;

% Initialize variables
% Note: Because of MATLAB's indexing, t=0 is index 1
r = zeros(num_nodes,3,sim_step+1);  % increments of del_t starting at 0
r(:,:,1) = r0;
v = zeros(num_nodes,3,sim_step);  % increments of del_t starting at del_t/2
KE = zeros(sim_step,1);
F_cable = zeros(num_nodes,3,sim_step);
F_rod = zeros(num_nodes,3,sim_step);
F_total = zeros(num_nodes,3,sim_step);
rod_length = zeros(num_rods,1,sim_step+1);

% Set mass and spring constants
m = 10000;      % mass per node
k_cable = 1;    % spring constant of the cables
k_rod = 1000;   % spring constant of the rods

% Set desired rest lengths
% rest_length = rand(num_cables,1).*original_cable_length;
rest_length = 0.9*original_cable_length;

%% Dynamic relaxation

for i = 1:sim_step

    % Find forces acting on each node
    F_cable(:,:,i) = findCableForce(r(:,:,i), cable_pair, num_nodes, ...
        k_cable, rest_length);
    [F_rod(:,:,i),rod_length(:,:,i)] = findRodForce(r(:,:,i), rod_pair, ...
        num_nodes, num_rods, k_rod, original_rod_length);
    F_total(:,:,i) = F_cable(:,:,i) + F_rod(:,:,i);
    
    % Find velocity at t + del_t/2, and velocity centered at t using
    % average of t+del_t/2 and t-del_t/2
    if i > 1
        v(:,:,i) = v(:,:,i-1) + del_t/m*F_total(:,:,i);
        v_center = (v(:,:,i) + v(:,:,i-1)) / 2;
    else
        % First velocity update is modified
        v(:,:,i) = del_t/m*F_total(:,:,i);
        v_center = zeros(size(v(:,:,1)));
    end
    
    % Find kinetic energy (KE)
    for j = 1:num_nodes
        KE(i) = KE(i) + 0.5*m*v_center(j,:)*v_center(j,:)';
    end
    
    % Kinetic damping: Reset velocity and KE to zero if peak was detected
    if i > 1
        if KE(i) - KE(i-1) < 0
            v(:,:,i) = 0;
            KE(i) = 0;
        end
    end
    
    % Update position
    r(:,:,i+1) = r(:,:,i) + v(:,:,i)*del_t;
    
end

% Store updated rod length
rod_length(:,:,end) = norm(r(rod_pair(1,1),:,end) - ...
    r(rod_pair(1,2),:,end));

%% Output results
% F_total
% rod_length

%% Plot results

% Kinetic energy
figure
plot(0:sim_step-1,KE);
xlabel('Time step')
ylabel('Kinetic energy')
grid on

% Tensegrity
figure
labels_on = 0;
plotTensegrity(r0, cable_pair, rod_pair, num_nodes, num_cables, ...
    num_rods, labels_on)
hold on
plotTensegrity(r(:,:,end), cable_pair, rod_pair, num_nodes, num_cables, ...
    num_rods, labels_on)
addForceToPlot(r(:,:,end),F_rod(:,:,end),'r')
addForceToPlot(r(:,:,end),F_cable(:,:,end),'m')