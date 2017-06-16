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

% System parameters
num_nodes = size(r0,1);
num_cables = size(cable_pairs,1);
num_rods = size(rod_pairs,1);

%% Plot tensegrity
% figure
% labels_on = 1;
% plotTensegrity(r0, cable_pairs, rod_pairs, labels_on)
% addCoordinateSystem(r0, rod_pairs)

%% Set up dynamic relaxation

% Set simulation variables
sim_steps = 10;
del_t = 0.01;

% Initialize position
r = zeros(num_nodes,3,sim_steps);
r(:,:,1) = r0;

% Initialize velocity and acceleration
% Note: Velocity (r_dot) will be recorded at every half step. Although not
% physically meaningful at t = 0, the point v(t-del_t/2) is the first entry
% to enable centered finite difference form of velocity to find
% acceleration. Then v(t), v(t+del_t/2), v(2*t), and so on.
r_dot = zeros(num_nodes,3,2*sim_steps);
r_dotdot = zeros(num_nodes,3,sim_steps);

% Initialize five state description of rods (x,y,z,phi,theta)
q(:,:,:) = zeros(num_rods,5,sim_steps);
q_dot(:,:,:) = zeros(size(q));

% Find first states
[q(:,:,1), q_dot(:,:,1)] = findState_rv(r(:,:,1), r_dot(:,:,1), rod_pairs);

%% Dynamic relaxation

% Set mass and spring constant
m = 1;      % mass per node
k = 10;     % spring constant

% Set desired rest lengths
% rest_lengths = rand(num_cables,1)
rest_lengths = 0.5*ones(num_cables,1);



% Initialize force variables
% F_cable = zeros(num_nodes,3,sim_steps);
% F_rod = zeros(num_nodes,3,sim_steps);

% for i = 1:sim_steps
% 
%     % Find cable force acting on each node
%     F_cable = findCableForce(r(:,:,i), cable_pairs, k, rest_lengths);
% 
%     % Find rod force acting on each node
%     F_rod = findRodForce(r(:,:,i), rod_pairs, q(:,:,i), q_dot(:,:,i), ...
%         q_dotdot(:,:,i), F_cable, m);
% 
%     % Find total force
%     F_total = F_cable + F_rod;
% 
%     % Find nodal positions, accelerations, and velocities
%     r_dotdot = F_total/m;
%     r_dot = r_dot_last + del_t/(2*m)*F_total;
%     r(:,:,i+1) = r(:,:,i) + r_dot*del_t;
%     
%     % Update states
%     q(:,:,i+1) = findState(r(:,:,i+1), rod_pairs);
%     q_dot(:,:,i+1) = findState(r_dot, rod_pairs);
%     q_dotdot(:,:,i+1) = findState(r_dotdot, rod_pairs);
% 
%     
%     % Find Kinetic Energy
% 
% end

%% Check forces
% F_cable
% F_rod
% 
% for i = 1:size(F_cable)
%     disp(' ')
%     disp(['Node ' num2str(i) ':'])
%     disp(['Cable force norm: ' num2str(norm(F_cable(i,:)))]);
%     disp(['rod force norm: ' num2str(norm(F_rod(i,:)))]);
% end