function intersect_found = checkRodIntersection(r, rod_pair, num_rods, ...
    rod_radius)
% This function checks if any two of the rods are intersected.
%
% The inputs are the following:
%   r: matrix of x,y,z, position of nodes
%   rod_pair: each row defines the node indices corresponding to that rod
%   num_rods: number of rods
%   rod_radius: radius of the rod
%
% The outputs are the following:
%   intersect_found: boolean indicating whether an intersection was found

% Find combinations of rods
comb = combnk(1:num_rods,2);

% Check for intersections
intersect_found = 0;
for i = 1:size(comb,1)
    
    % Find distance of closest point to both rods
    rod_A = comb(i,1);
    rod_B = comb(i,2);
    node_A_start = rod_pair(rod_A,1);
    node_A_end = rod_pair(rod_A,2);    
    node_B_start = rod_pair(rod_B,1);
    node_B_end = rod_pair(rod_B,2);
    PA = [r(node_A_start,:); r(node_B_start,:)];
    PB = [r(node_A_end,:); r(node_B_end,:)];
    [~,distance] = lineIntersect3D(PA,PB);
    
    % Check if the distance is less than the rod's radius
    for j = 1:length(distance)
        if distance(j) < rod_radius
            intersect_found = 1;
        end
    end
    
end