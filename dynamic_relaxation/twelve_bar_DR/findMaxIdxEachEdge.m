function max_idx = findMaxIdxEachEdge(results)
% This function finds the indices of the cable combinations that resulted
% in the greatest distance of the COG's escape from the supporting polygon
% relative to each edge. Designed to be used with the combinatoric approach
% script.
%
% Inputs:
%   results: structure of results with fields edge_distance, and 
%       edge_closest
%   num_act_cables: the number of actuated cables
%
% Outputs:
%   max_idx: each entry is the index corresponding to the maximum distance
%       of escape of the COG from each edge of the polygon (entry is 0 if
%       no escape from a given edge)

% Grab number of edges
num_edges = size(results.ground_face,2);

% Initialize output
max_idx = zeros(num_edges,1);

% Fill in results summary
for i = 1:num_edges
    
    % Find results that moved the COG past the current edge
    idx = find(results.edge_closest == i);
    if ~isempty(idx)
        
        % Get information for maximum displacement off this edge
        [~,max_local_idx] = max(results.edge_distance(idx));
        max_idx(i) = idx(max_local_idx);
%         edge_cable_combo = results.cable_combos(max_idx(i));
%         intersection = results.intersect(max_idx(i));
        
%         % Print out information
%         fprintf(['\nThe COG escapes from this edge by a distance of ' ...
%             '%.2f using the following actuated cable(s):\n'], max_distance)
%         disp(edge_cable_combo)
%         if intersection == 1
%             fprintf('\nHowever, an intersection was found on this edge.\n')
%         end
        
    end
end




% % Plotting
% % Final tensegrity
%     figure
%     plotTensegrity(r_rot, cable_pair, rod_pair, labels_on, color)
%     % addForceToPlot(r(:,:,end),F_total(:,:,end),'g')  % plot total forces
%     hold on
%     r_rot_ground = r_rot(ground_face,:);
%     scatter3(r_rot_ground(:,1),r_rot_ground(:,2),r_rot_ground(:,3),'Filled','b')
%     hold on
%     scatter3(COG(1),COG(2),COG(3),'Filled',color)
%     hold on
%     plot3([COG(1); COG(1)], [COG(2); COG(2)], [COG(3); ...
%         min(r_rot(:,3))],['--' color])
%     title(['Cables used: ' num2str(num_comb) ', Combination #: ' num2str(row(i))])