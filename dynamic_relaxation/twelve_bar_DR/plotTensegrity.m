function plotTensegrity(r, cable_pair, rod_pair, labels_on, cable_format)
% Generic tensegrity plotter.
%
% The inputs are the following:
%   r: matrix of x,y,z, position of nodes
%   cable_pair: each row defines the node indices corresponding to that
%     cable
%   rod_pair: each row defines the node indices corresponding to that bar
%   labels_on: boolean to add labels of bars and cables to plot
%   cable_format: character string made from elements of plotting columns
%       (e.g. '-k')

% Grab number of nodes, cables, and rods
num_nodes = size(r,1);
num_cables = size(cable_pair,1);
num_rods = size(rod_pair,1);

% Plot nodes
for i = 1:num_nodes
    hold on
    scatter3(r(i,1),r(i,2),r(i,3),'k','filled')
    if labels_on == 1
        text(r(i,1),r(i,2),r(i,3),...
            ['  ' num2str(i-1) '  '],'FontWeight','Bold')
    end
end 
axis equal
xlabel('x (cm)')
ylabel('y (cm)')
zlabel('z (cm)')
view(45,15)
% set(gca,'FontSize',16)
grid on

% Plot cables
for i = 1:num_cables
    hold on
    cable_x = [r(cable_pair(i,1),1); r(cable_pair(i,2),1)];
    cable_y = [r(cable_pair(i,1),2); r(cable_pair(i,2),2)];
    cable_z = [r(cable_pair(i,1),3); r(cable_pair(i,2),3)];
    plot3(cable_x,cable_y,cable_z,'b','LineWidth',1.5,'Color',cable_format)
    if labels_on == 1
        text(sum(cable_x)/2,sum(cable_y)/2,sum(cable_z)/2,...
            ['  ' num2str(i-1) '  '],'Color',cable_format)
    end
end

% Plot bars
for i = 1:num_rods
    hold on
    bar_x = [r(rod_pair(i,1),1); r(rod_pair(i,2),1)];
    bar_y = [r(rod_pair(i,1),2); r(rod_pair(i,2),2)];
    bar_z = [r(rod_pair(i,1),3); r(rod_pair(i,2),3)];
    plot3(bar_x,bar_y,bar_z,'Color',[0.7 0.7 0.7],'LineWidth',2.5)
%     plot3(bar_x,bar_y,bar_z,'k','LineWidth',2.5)
    if labels_on == 1
        text(sum(bar_x)/2,sum(bar_y)/2,sum(bar_z)/2,...
            ['  ' num2str(i+9-1) '  '])
    end
end