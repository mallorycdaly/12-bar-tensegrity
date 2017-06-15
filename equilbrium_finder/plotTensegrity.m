function plotTensegrity(nodes, cable_pairs, bar_pairs, labels_on)
% Generic tensegrity plotter. The inputs are the following:
%   nodes: matrix of x,y,z, position of nodes
%   cable_pairs: each row defines the node indices corresponding to that
%     cable
%   bar_pairs: each row defines the node indices corresponding to that bar
%   labels_on: boolean to add labels of bars and cables to plot

% Plot nodes
scatter3(nodes(:,1),nodes(:,2),nodes(:,3),'k')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

% Grab number of cables and bars
num_cables = size(cable_pairs,1);
num_bars = size(bar_pairs,1);

% Plot cables
for i = 1:num_cables
    hold on
    cable_x = [nodes(cable_pairs(i,1),1); nodes(cable_pairs(i,2),1)];
    cable_y = [nodes(cable_pairs(i,1),2); nodes(cable_pairs(i,2),2)];
    cable_z = [nodes(cable_pairs(i,1),3); nodes(cable_pairs(i,2),3)];
    plot3(cable_x,cable_y,cable_z,'k','LineWidth',1.5)
    if labels_on == 1
        text(sum(cable_x)/2,sum(cable_y)/2,sum(cable_z)/2,['  ' num2str(i) '  '])
    end
end

% Plot bars
for i = 1:num_bars
    hold on
    bar_x = [nodes(bar_pairs(i,1),1); nodes(bar_pairs(i,2),1)];
    bar_y = [nodes(bar_pairs(i,1),2); nodes(bar_pairs(i,2),2)];
    bar_z = [nodes(bar_pairs(i,1),3); nodes(bar_pairs(i,2),3)];
    plot3(bar_x,bar_y,bar_z,'k','LineWidth',2.5)
    if labels_on == 1
        text(sum(bar_x)/2,sum(bar_y)/2,sum(bar_z)/2,['  ' num2str(i+9) '  '])
    end
end