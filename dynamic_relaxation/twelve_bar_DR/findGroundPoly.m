function ground_poly = findGroundPoly(r, ground_face)
% This function compiles the nodes of the ground polygon based on the
% ground face.
%
% The inputs are the following:
%   r: matrix of x,y,z positions of the nodes for the given face
%   ground_face: vector of node indices for the given face
%
% The outputs are the following:
%   ground_poly: x,y,z positions of nodes that form the ground polygon

ground_poly = zeros(length(ground_face),3);
for i = 1:length(ground_face)
    ground_poly(i,:) = r(ground_face(i),:);
end