function [r, v, KE, F_cable, F_rod, F_total, L_rod] = ...
    dynamicRelaxation(r0, cable_pair, rod_pair, m, k_spring, L0_spring, ...
    k_rod, L0_rod, rest_lengths, sim_steps, del_t)
% Runs dynamic relaxation from intial position to final configuration,
% based on number of simulation steps.
%
% The inputs are the following:
%   r0: matrix of initial x,y,z, position of nodes
%   cable_pair: each row defines the node indices corresponding to that
%       cable
%   rod_pair: each row defines the node indices corresponding to that rod
%   rod_radius: radius of the rod
%   m: mass of each node
%   k_spring: spring constant of the springs
%   L0_spring: initial length of the springs
%   k_rod: spring constant of the rods
%   L0_rod: initial length of the rods
%   rest_lengths: lengths of strings in the cable, defined as the cable 
%       rest lengths, such that spring length can be found as total cable 
%       length (node-to-node distance) minus rest length
%   sim_steps: number of steps in the simulation
%   del_t: time incrementation
%
% The outputs are the following:
%   r: 3D array of x,y,z, position of nodes across simulation steps
%   v: 3D array of x,y,z, velocity of nodes across simulation steps
%   KE: vector of kinetic energy across simulation steps
%   F_cable: 3D array of cable force vectors across simulation steps
%   F_rod: 3D array of rod force vectors across simulation steps
%   F_total: 3D array of total force vectors across simulation steps
%   L_rod: 3D array of rod lengths across simulation steps

% Grab number of nodes and rods
num_nodes = size(r0,1);
num_rods = size(rod_pair,1);

% Initialize variables
% Note: All variables except velocity are incremented by del_t starting
% from t = 0 (index 1 in MATLAB). Velocity is incremented by del_t starting
% from t = -del_t/2 (index 1).
r = zeros(num_nodes,3,sim_steps+1);      % nodal positions
r(:,:,1) = r0;                           % initial nodal positions known
v = zeros(num_nodes,3,sim_steps+1);      % nodal velocities
KE = zeros(sim_steps,1);                 % kinetic energy
F_cable = zeros(num_nodes,3,sim_steps);  % cable force at nodes
F_rod = zeros(num_nodes,3,sim_steps);    % rod force at nodes
F_total = zeros(num_nodes,3,sim_steps);  % total force at nodes
L_rod = zeros(num_rods,1,sim_steps+1);   % rod length

% Run dynamic relaxation
for i = 1:sim_steps

    % Find forces acting on each node
    F_cable(:,:,i) = findCableForce(r(:,:,i), cable_pair, num_nodes, ...
        k_spring, L0_spring, rest_lengths);
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

% % Check that forces converged to near zero and throw warning if not
% F_total_rounded = round(F_total(:,:,end),6);
% if any(any(F_total_rounded ~= 0))
%     warning(['Total rod forces did not converge to zero with a ' ...
%         'precision of 6 decimal places'])
% end