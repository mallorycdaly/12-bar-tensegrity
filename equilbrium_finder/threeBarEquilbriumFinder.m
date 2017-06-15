%% Test optimization-based equilbrium finder using three strut prism
clear; close all

%% Set up three-strut prism tensegrity
% Using prism numbering scheme from Fig. 10 of Review of Form-Finding
% Methods for Tensegrity Structures by A. Tibert.

% Arbitrarily assign design parameters
edge_length = 1;
height = 1;
rotation_angle = -25;  % degrees

% Find nodes based on design parameters
[nodes, cable_pairs, bar_pairs] = formThreeBar(edge_length, height, ...
    rotation_angle);

%% Plot tensegrity
figure
labels_on = 1;
plotTensegrity(nodes, cable_pairs, bar_pairs, labels_on)
addCoordinateSystem(nodes, bar_pairs)

%% Find initial states
q0 = findState(nodes, bar_pairs)

%% Dynamic relaxation


































%%
% % Following M. Barnes' thesis: Form finding and analysis of tension space
% % structures by dynamic relaxation
% %
% % Governing equation: F_xi(t) = m_i*v_dot_xi(t) + c_i*v_xi(t) = k*x
% %   F_xi = residual force at node i in direction x at time t
% %   m_i = mass at node i
% %   c_i = damping constant at node i
% %   v_xi, V_dot_xi = velocity and acceleration of node i in direction x at
% %       time t
% 
% % Define mass and damping constant per node
% m = 1;
% c = 0.5;
% del_t = 0.01;
% A = 1/(m/del_t + c/2);
% B = (m/del_t - c/2)/(m/del_t + c/2);
% k = 10;
% 
% % Set rest lengths
% rest_lengths = 0.5*ones(num_cables,1);
% 
% % Initialize V and V_dot
% V = zeros(size(nodes));
% V_dot = zeros(size(nodes));
% 
% % for i = 1:numNodes
% 
% % Find initial forces

    