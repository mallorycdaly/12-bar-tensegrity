%% Test optimization-based equilbrium finder using three strut prism
% Using prism numbering scheme from Fig. 10 of Review of Form-Finding
% Methods for Tensegrity Structures by A. Tibert.
clear all
close all

%% Set up tensegrity

% System definition
numCables = 9;
numBars = 3;

% Arbitrarily assign design parameters
edgeLength = 1;
height = 1;

% Base triangle centered at (0,0,0)
node4 = [edgeLength/2 -sqrt(3)/6*edgeLength 0];
node5 = [0 sqrt(3)/3*edgeLength 0];
node6 = [-edgeLength/2 -sqrt(3)/6*edgeLength 0];

% Rotation angle in degrees
rotAngle = -45;
rotM = [cosd(rotAngle) -sind(rotAngle) 0; sind(rotAngle) cosd(rotAngle) 0; 0 0 1];

% Rotate base triangle to form top triangle
node1 = node4*rotM + [0 0 height];
node2 = node5*rotM + [0 0 height];
node3 = node6*rotM + [0 0 height];

% Assemble nodes in matrix
Nodes = [node1; node2; node3; node4; node5; node6];

% Cables and bar pairs
cablePairs = [1 2; 2 3; 1 3; 1 4; 2 5; 3 6; 4 5; 5 6; 4 6];
barPairs = [1 6; 2 4; 3 5];

%% Plot tensegrity
figure

% Plot nodes
scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3))
axis equal

% Plot cables
for i = 1:numCables
    hold on
    cable_x = [Nodes(cablePairs(i,1),1); Nodes(cablePairs(i,2),1)];
    cable_y = [Nodes(cablePairs(i,1),2); Nodes(cablePairs(i,2),2)];
    cable_z = [Nodes(cablePairs(i,1),3); Nodes(cablePairs(i,2),3)];
    plot3(cable_x,cable_y,cable_z,'Color',[0.7 0.7 0.7],'LineWidth',1.5)
%     text(sum(cable_x)/2,sum(cable_y)/2,sum(cable_z)/2,['  ' num2str(i) '  '])
end

% Plot bars
for i = 1:numBars
    hold on
    bar_x = [Nodes(barPairs(i,1),1); Nodes(barPairs(i,2),1)];
    bar_y = [Nodes(barPairs(i,1),2); Nodes(barPairs(i,2),2)];
    bar_z = [Nodes(barPairs(i,1),3); Nodes(barPairs(i,2),3)];
    plot3(bar_x,bar_y,bar_z,'Color',[0.7 0.7 0.7],'LineWidth',2.5)
%     text(sum(bar_x)/2,sum(bar_y)/2,sum(bar_z)/2,['  ' num2str(i+9) '  '],'Color','red')
end

%% 