%% 12-Bar Tensegrity Cube
% From SW model 'C:\Users\Mallory\Google Drive\Tensegrity Research at NASA
% Ames and UC Berkeley\12-Bar Tensegrity\12-Bar Structures\Octahedron\
% octahedron adjustable dims.sldprt'
% See also 'C:\Users\Mallory\Google Drive\Tensegrity Research at NASA Ames
% and UC Berkeley\12-Bar Tensegrity\12-Bar Structures\Node Finding\node
% numbers.ppx'
% Author: Mallory Daly
% Created: 2/22/17
% Last Updated: 4/3/17
close all

%% Scaling Factor
% Set scaling factor to 0.1 to mimic NTRT conversion to dm
scalingFactor = 0.1;

%% Node Positions
% Rod length = 45 cm

% From SW model
% n0  = [-33.37 5.83 5.82];
% n1  = [-11.63 -31.81 -5.82];
% n2  = [-19.49 25.98 -11.25];
% n3  = [19.49 25.98 11.25];
% n4  = [-2.53 30.36 -15.91];
% n5  = [19.2 -7.28 -27.56];
% n6  = [-3.3 20.27 -27.56];
% n7  = [-25.03 -17.37 -15.91];
% n8  = [22.5 12.99 -22.5];
% n9  = [22.5 12.99 22.5];
% n10 = [30.74 -12.99 -8.24];
% n11 = [-8.24 -12.99 -30.74];
% n12 = [25.03 -17.37 15.91];
% n13 = [3.30 20.27 27.56];
% n14 = [11.63 -31.81 5.82];
% n15 = [33.37 5.83 -5.82];
% n16 = [0 -25.98 22.5];
% n17 = [0 -25.98 -22.5];
% n18 = [-19.2 -7.28 27.56];
% n19 = [2.53 30.36 15.91];
% n20 = [-22.5 12.99 22.5];
% n21 = [-22.5 12.99 -22.5];
% n22 = [-30.74 -12.99 8.24];
% n23 = [8.24 -12.99 30.74];

% Measured by hand
n0 = [0 0 0];
n1 = [5.4 3.5 18.4];
n2 = [23.2 0 13.8];
n3 = [34.9 3.5 0];
n4 = [23.2 0 -13.8];
n5 = [5.4 3.5 -18.4];
n6 = [-6.2 18 -18.1];
n7 = [-12.2 14 0];
n8 = [9.7 18 27];
n9 = [27.2 14 23.8];
n10 = [41.8 18 -8.5];
n11 = [28.9 14 -22.5];
n12 = [6.8 27 -24.6];
n13 = [-12.4 31 2.5];
n14 = [-4 27 18.2];
n15 = [29.5 31 21.3];
n16 = [39.5 27 4.6];
n17 = [26 31 -24.3];
n18 = [0 41.5 -15.5];
n19 = [0 45 0];
n20 = [3.6 41.5 20.1];
n21 = [23.2 45 13.8];
n22 = [36.5 41.5 0];
n23 = [23.2 45 -13.8];

% Scale nodes by compiling as matrix and scaling matrix
nodeRod = scalingFactor*[n0; n1; n2; n3; n4; n5; n6; n7; n8; n9; n10; n11; 
    n12; n13; n14; n15; n16; n17; n18; n19; n20; n21; n22; n23];

% Plot nodes
figure
% scatter3(nodeRod(:,1),nodeRod(:,2),nodeRod(:,3),'filled')
title('12-Bar Tensegrity Octahedron')
% viewlim = [-3 3];
% xlim(viewlim)
% ylim(viewlim)
% zlim(viewlim)
grid on
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

%% Bar Positions
% Rod length = 45 cm

% Numbering convention in 'node numbers.ppx'
% b0  = [n0;  n1];
% b1  = [n2;  n3];
% b2  = [n4;  n5];
% b3  = [n6;  n7];
% b4  = [n8;  n9];
% b5  = [n10; n11];
% b6  = [n12; n13];
% b7  = [n14; n15];
% b8  = [n16; n17];
% b9  = [n18; n19];
% b10 = [n20; n21];
% b11 = [n22; n23];

% Numbering convention written on octahedron
b0  = [n0;  n19];
b1  = [n2;  n21];
b2  = [n4;  n23];
b3  = [n6;  n3];
b4  = [n13;  n9];
b5  = [n18; n16];
b6  = [n1; n10];
b7  = [n7; n17];
b8  = [n14; n22];
b9  = [n20; n12];
b10 = [n15; n11];
b11 = [n8; n5];

% Scale bars by compiling as matrix and scaling matrix
pairRod = scalingFactor*[b0; b1; b2; b3; b4; b5; b6; b7; b8; b9; b10; b11];

% Add bars to plot
colors = 'bgrcmy';
for i = 1:2:23
    hold on
    plot3(pairRod(i:i+1,1),pairRod(i:i+1,2),pairRod(i:i+1,3),'Color',[0.7 0.7 0.7],'LineWidth',3.5)
    grid on
end
xlabel('x (dm)')
ylabel('y (dm)')
zlabel('z (dm)')

%% Cable Positions
% Cables were not designed to be equal length in SW model

% Numbering convention in 'node numbers.ppx'
% c0 = [n0; n2];
% c1 = [n0; n20];
% c2 = [n0; n22];
% c3 = [n2; n4];
% c4 = [n2; n21];
% c5 = [n4; n6];
% c6 = [n4; n19];
% c7 = [n6; n8];
% c8 = [n6; n21];
% c9 = [n8; n5];
% c10 = [n8; n15];
% c11 = [n5; n10];
% c12 = [n5; n17];
% c13 = [n10; n12];
% c14 = [n10; n15];
% c15 = [n12; n14]; 
% c16 = [n12; n23];
% c17 = [n14; n1];
% c18 = [n14; n16];
% c19 = [n16; n18];
% c20 = [n16; n23];
% c21 = [n18; n20];
% c22 = [n18; n22];
% c23 = [n20; n13];
% c24 = [n21; n11];
% c25 = [n15; n3];
% c26 = [n17; n11];
% c27 = [n17; n1];
% c28 = [n1; n7];
% c29 = [n22; n7];
% c30 = [n19; n3];
% c31 = [n19; n13];
% c32 = [n11; n7];
% c33 = [n9; n3];
% c34 = [n9; n23];
% c35 = [n9; n13];

% Numbering convention written on octahedron
c0 = [n0; n1];
c1 = [n0; n5];
c2 = [n0; n7];
c3 = [n1; n2];
c4 = [n1; n8];
c5 = [n2; n3];
c6 = [n2; n9];
c7 = [n3; n4];
c8 = [n3; n10];
c9 = [n4; n5];
c10 = [n4; n11];
c11 = [n5; n6];
c12 = [n6; n7];
c13 = [n6; n12];
c14 = [n7; n13];
c15 = [n8; n9]; 
c16 = [n8; n14];
c17 = [n9; n15];
c18 = [n10; n11];
c19 = [n10; n16];
c20 = [n11; n17];
c21 = [n12; n17];
c22 = [n12; n18];
c23 = [n13; n14];
c24 = [n13; n19];
c25 = [n14; n20];
c26 = [n15; n16];
c27 = [n15; n21];
c28 = [n16; n22];
c29 = [n17; n23];
c30 = [n18; n19];
c31 = [n18; n23];
c32 = [n19; n20];
c33 = [n20; n21];
c34 = [n21; n22];
c35 = [n22; n23];

% Scale cables by compiling as matrix and scaling matrix
pairStructureCable = scalingFactor*[c0; c1; c2; c3; c4; c5; c6; c7; c8; c9; c10;
    c11; c12; c13; c14; c15; c16; c17; c18; c19; c20; c21; c22; c23; c24;
    c25; c26; c27; c28; c29; c30; c31; c32; c33; c34; c35];

% Add cables to plot
for i = 1:2:71
    hold on
    plot3(pairStructureCable(i:i+1,1),pairStructureCable(i:i+1,2),pairStructureCable(i:i+1,3),'r','LineWidth',1.25)
end

%% Add payload to center

% Find center of tensegrity
center = round(mean(nodeRod),4);

% Define edge length of payload
edgeLength = 0.8; %0.8 dm = 80 mm

% Form and plot cubical payload
x = ([0 1 1 0 0 0; 1 1 0 0 1 1; 1 1 0 0 1 1; 0 1 1 0 0 0] - 0.5)*edgeLength + center(1);
y = ([0 0 1 1 0 0; 0 1 1 0 0 0; 0 1 1 0 1 1; 0 0 1 1 1 1] - 0.5)*edgeLength + center(2);
z = ([0 0 0 0 0 1; 0 0 0 0 0 1; 1 1 1 1 0 1; 1 1 1 1 0 1] - 0.5)*edgeLength + center(3);
for i = 1:6
    v = [x(:,i),y(:,i),z(:,i)];
    f = [1 2 3 4];
    hold on
    h = patch('Faces',f,'Vertices',v);
    set(h,'edgecolor',[0.7 0.7 0.7],'facecolor',[0.85 0.85 0.85])
end

%% Find and plot nodes for center of each face of payload

% Define nodes
nodePayload = [mean(x); mean(y); mean(z)]';

% Add to plot1
hold on
scatter3(mean(x),mean(y),mean(z),'k.')

%% Find nodes for center of each rod

nodeRodCenter = zeros(12,3);
for i = 1:2:23
    j = (i+1)/2;
    nodeRodCenter(j,:) = mean(pairRod(i:i+1,:));
end

%% Find and add payload suspension cables

% Find three closest rod centers to each payload face center
pairPayloadCable = zeros(2*3*6,3);
for i = 1:size(nodePayload,1)
    diff = nodePayload(i,:)-nodeRodCenter;
    dist = sqrt(sum(diff.^2,2));
    [~, index] = sort(dist);
    pairPayloadCable(6*(i-1)+1:6*i,:) = ...
        [nodePayload(i,:); nodeRodCenter(index(1),:);
         nodePayload(i,:); nodeRodCenter(index(2),:);
         nodePayload(i,:); nodeRodCenter(index(3),:)];
end

% Add cables to plot
for i = 1:2:(size(pairPayloadCable,1)-1)
    hold on
    plot3(pairPayloadCable(i:i+1,1),pairPayloadCable(i:i+1,2),pairPayloadCable(i:i+1,3),'r','LineWidth',1)
end