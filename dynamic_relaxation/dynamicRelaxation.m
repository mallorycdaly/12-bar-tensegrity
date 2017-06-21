function [r_new, v_new, KE, F_cable, F_rod, F_total, L_rod] = ...
    dynamicRelaxation(r, v, rod_pair, cable_pair, m, k_rod, k_spring, L0_spring, L0_rod, rest_length, sim_idx, del_t)

% Grab number of nodes and rods
num_nodes = size(r,1);
num_rods = size(rod_pair,1);

% Find forces acting on each node
F_cable = findCableForce(r, cable_pair, num_nodes, k_spring, L0_spring, rest_length);
[F_rod, L_rod] = findRodForce(r, rod_pair, num_nodes, num_rods, k_rod, L0_rod);
F_total = F_cable + F_rod;

% Find velocity at t + del_t/2, and velocity centered at t using
% average of t+del_t/2 and t-del_t/2

    % First velocity update is modified
    if sim_idx == 1
        v_new = F_total/(m/del_t);            
        v_center = zeros(size(v));
    else
        v_new = v + F_total/(m/del_t);
        v_center = (v + v_new) / 2;
    end

% Find kinetic energy (KE)
KE = 0;
for j = 1:num_nodes
    KE = KE + 0.5*m*v_center(j,:)*v_center(j,:)';
end

% Kinetic damping: Reset velocity and KE to zero if peak was detected
% and flag for form finding process to restart
if sim_idx > 1
    if KE(sim_idx) - KE(sim_idx-1) < 0
        % Reset velocity and KE
        v_new(:,:,sim_idx) = 0;
        KE(sim_idx) = 0;
    end
end

% Update position
r_new = r + v_new*del_t;
