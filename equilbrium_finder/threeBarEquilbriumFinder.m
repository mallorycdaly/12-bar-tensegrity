%% Test optimization-based equilbrium finder using three strut prism
clear; close all

%% Set up three-strut prism tensegrity
% Using prism numbering scheme from Fig. 10 of Review of Form-Finding
% Methods for Tensegrity Structures by A. Tibert.

% Arbitrarily assign design parameters
edge_length = 1;
height = 1;
rotation_angle = -25;  % degrees

% Find initial position of nodes (r0) and indices corresponding to cables
% (cable_pairs) and rods (rod_pairs)
[r0, cable_pairs, rod_pairs] = formThreeBar(edge_length, height, ...
    rotation_angle);

% Original rod length, before any deformation
original_rod_length = norm(r0(rod_pairs(1,1),:) - r0(rod_pairs(1,2),:));

% System parameters
num_nodes = size(r0,1);
num_cables = size(cable_pairs,1);
num_rods = size(rod_pairs,1);

%% Plot original tensegrity
% figure
% labels_on = 0;
% plotTensegrity(r0, cable_pairs, rod_pairs, labels_on)
% addCoordinateSystem(r0, rod_pairs)

%% Set up dynamic relaxation

% Set simulation variables
sim_steps = 1000;
del_t = 0.0001;

% Initialize position, stored at increments of t+del_t
r = zeros(num_nodes,3,sim_steps+1);
r(:,:,1) = r0;

% Initialize velocity, stored at increments of t+del_t/2
v = zeros(num_nodes,3,sim_steps+1);

% Initialize kinetic energy
KE = zeros(sim_steps+1,1);

% Set mass and spring constants
m = 1;          % mass per node
k_cable = 10;   % spring constant of the cables
k_rod = 10000;  % spring constant of the rods

% Set desired rest lengths
rest_lengths = rand(num_cables,1)
% rest_lengths = 0.5*ones(num_cables,1);

% % Initialize five state description of rods (x,y,z,phi,theta)
% q(:,:,:) = zeros(num_rods,5,sim_steps);
% q_dot(:,:,:) = zeros(size(q));
% 
% % Find first states
% [q(:,:,1), q_dot(:,:,1)] = findState_rv(r(:,:,1), v(:,:,1), rod_pairs);

%% Dynamic relaxation

for i = 1:sim_steps

    % Find forces acting on each node
    F_cable = findCableForce(r(:,:,i), cable_pairs, k_cable, rest_lengths);
    F_rod = findRodForce(r(:,:,i), rod_pairs, k_rod, original_rod_length);
    F_total = F_cable + F_rod;

    % Update velocity
    v(:,:,i+1) = v(:,:,i) + del_t/m*F_total;

    % Update position
    r(:,:,i+1) = r(:,:,i) + v(:,:,i+1)*del_t;

    % Find kinetic energy, KE = 1/2*m*v^2
    v_center = (v(:,:,i)+v(:,:,i+1))/2;
    for j = 1:num_nodes
        KE(i+1) = KE(i+1) + 0.5*m*v_center(j,:)*v_center(j,:)';
    end
    
    % Kinetic damping: reset velocities to zero if peak was detected
    if KE(i+1) - KE(i) < 0
        v(:,:,i+1) = 0;
    end
    
end

%% Plot results

% Kinetic energy
figure
plot(0:sim_steps,KE);
xlabel('Time step')
ylabel('Kinetic energy')
grid on

% Tensegrity
figure
labels_on = 0;
plotTensegrity(r(:,:,end), cable_pairs, rod_pairs, labels_on)