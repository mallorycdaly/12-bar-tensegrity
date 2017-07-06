function [normal_vec, COG_plane] = findAvgNormalVector(r, ground_face)
% This function finds the unit average normal vector to a polygon that is
% expected to be non-planar.
%
% The inputs are the following:
%   r: matrix of x,y,z positions of the nodes for the given face
%   ground_face: vector of node indices for the given face
%
% The outputs are the following:
%   normal_vec = the average normal vector found from the nodes of the
%       ground face

face_nodes = r(ground_face,:);
COG_plane = mean(face_nodes,1);

wrapN = @(x,N)(1+mod(x,N));

indiv_normal_vec = zeros(length(ground_face),3);
for i = 1:length(ground_face)
    nodeA_idx = wrapN(i,length(ground_face));
    nodeB_idx = wrapN(i+1,length(ground_face));
    nodeC_idx = wrapN(i+2,length(ground_face));
    nodeA = r(ground_face(nodeA_idx),:);
    nodeB = r(ground_face(nodeB_idx),:);
    nodeC = r(ground_face(nodeC_idx),:);
    indiv_normal_vec(i,:) = cross(nodeA-nodeB,nodeC-nodeB);
end

normal_vec = mean(indiv_normal_vec,1)/norm(mean(indiv_normal_vec,1));