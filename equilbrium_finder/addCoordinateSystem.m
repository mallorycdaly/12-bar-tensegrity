function addCoordinateSystem(nodes, bar_pairs)
% This function adds the world Cartesian coordinate system to a plot, as
% well as the local Cartesian coordinate system on each rod based on the
% translation:
%   r = x*E_1 + y*E_2 + z*E_3
% 
% The local coordinate system is based on the Fig. 4 of Kyunam Kim's
% "Robust Learning of Tensegrity Robot Control for Locomotion through
% Form-Finding."
% 
% The inputs are the following:
%   nodes: matrix of x,y,z, position of nodes
%   bar_pairs: each row defines the node indices corresponding to that bar

% Grab number of bars
num_bars = size(bar_pairs,1);

% World basis, Cartesian
I = 0.2*eye(3);
str = 'xyz';
for i = 1:3
    hold on
    plot3([0 I(i,1)],[0 I(i,2)],[0 I(i,3)],'r','LineWidth',1.5)
    text(I(i,1),I(i,2),I(i,3),[' ' str(i) ' '],'Color','red')
end

% Reorder nodes by bar pairs and so that higher (z) point is second
% Odd nodes: lower value of z
% Even nodes: higher value of z
nodes_bar = zeros(size(nodes));
for i = 1:num_bars
    if (nodes(bar_pairs(i,1),3) > nodes(bar_pairs(i,2),3))
        nodes_bar(2*i-1:2*i,:) = [nodes(bar_pairs(i,2),:)
            nodes(bar_pairs(i,1),:)];
    else
        nodes_bar(2*i-1:2*i,:) = [nodes(bar_pairs(i,1),:)
            nodes(bar_pairs(i,2),:)];
    end
end

% Plot local coordinate systems
for i = 1:3
    r = nodes_bar(2*i-1,:);
    for j = 1:3
        hold on
        plot3([0 I(j,1)]+r(1),[0 I(j,2)]+r(2),[0 I(j,3)]+r(3),'r','LineWidth',1)
        text(I(j,1)+r(1),I(j,2)+r(2),I(j,3)+r(3),[' ' str(j) ' '],'Color','red')
    end
end