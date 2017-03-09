%% 12-Bar Tensegrity Cube
% From SW model 'C:\Users\Mallory\Google Drive\Tensegrity Research at NASA
% Ames and UC Berkeley\12-Bar Tensegrity\12-Bar Structures\Cube\cube
% adjustable dims.sldprt'
% See also 'C:\Users\Mallory\Google Drive\Tensegrity Research at NASA Ames
% and UC Berkeley\12-Bar Tensegrity\12-Bar Structures\Node Finding\node
% numbers.ppx'
% Author: Mallory Daly
% Created: 2/22/17
% Last Updated: 2/22/17
close all

%% Scaling Factor
% Set scaling factor to 0.1 to mimic NTRT conversion to dm
scalingFactor = 0.1;

%% Node Positions
% Nodes created by starting from front view and moving clockwise from
% leftmost top node, and pairing nodes with bars
% Rod length = 45 cm
n0  = [-20.79 20.79 8.61];
n1  = [-8.61 -20.79 20.79];
n2  = [-20.79 20.79 -8.61];
n3  = [20.79 8.61 -20.79];
n4  = [-8.61 20.79 -20.79];
n5  = [-20.79 -20.79 -8.61];
n6  = [8.61 20.79 -20.79];
n7  = [20.79 8.61 20.79];
n8  = [20.79 20.79 -8.61];
n9  = [8.61 -20.79 -20.79];
n10 = [20.79 20.79 8.61];
n11 = [-20.79 8.61 20.79];
n12 = [8.61 20.79 20.79];
n13 = [20.79 -20.79 8.61];
n14 = [-8.61 20.79 20.79];
n15 = [-20.79 8.61 -20.79];
n16 = [-20.79 -20.79 8.61];
n17 = [20.79 -8.61 20.79];
n18 = [-8.61 -20.79 -20.79];
n19 = [-20.79 -8.61 20.79];
n20 = [20.79 -20.79 -8.61];
n21 = [-20.79 -8.61 -20.79];
n22 = [8.61 -20.79 20.79];
n23 = [20.79 -8.61 -20.79];

% Scale nodes by compiling as matrix and scaling matrix
nMatrix = scalingFactor*[n0; n1; n2; n3; n4; n5; n6; n7; n8; n9; n10; n11; 
    n12; n13; n14; n15; n16; n17; n18; n19; n20; n21; n22; n23];

% Plot nodes
figure
scatter3(nMatrix(:,1),nMatrix(:,2),nMatrix(:,3),'filled')
title('12-Bar Tensegrity Cube')
axis equal

%% Bar Positions
% Rod length = 45 cm
b0  = [n0;  n1];
b1  = [n2;  n3];
b2  = [n4;  n5];
b3  = [n6;  n7];
b4  = [n8;  n9];
b5  = [n10; n11];
b6  = [n12; n13];
b7  = [n14; n15];
b8  = [n16; n17];
b9  = [n18; n19];
b10 = [n20; n21];
b11 = [n22; n23];

% Scale bars by compiling as matrix and scaling matrix
barMatrix = scalingFactor*[b0; b1; b2; b3; b4; b5; b6; b7; b8; b9; b10; b11];

% Add bars to plot
for i = 1:2:23
    hold on
    plot3(barMatrix(i:i+1,1),barMatrix(i:i+1,2),barMatrix(i:i+1,3),'LineWidth',3.5)
end
xlabel('x (dm)')
ylabel('y (dm)')
zlabel('z (dm)')

%% Cable Positions
% Cables designed in SW model to be approximately equal length
c0	= [n0; n2];
c1	= [n0; n14];
c2	= [n0; n11];
c3	= [n2; n4];
c4	= [n2; n15];
c5	= [n4; n15];
c6	= [n4; n6];
c7	= [n6; n8];
c8	= [n6; n3];
c9	= [n8; n3];
c10	= [n3; n23];
c11	= [n8; n10];
c12	= [n10; n12];
c13	= [n10; n7];
c14	= [n7; n12];
c15	= [n14; n12];
c16	= [n14; n11];
c17	= [n23; n20];
c18	= [n9; n23];
c19	= [n7; n17];
c20	= [n17; n13];
c21	= [n13; n22];
c22	= [n22; n17];
c23	= [n22; n1];
c24	= [n1; n16];
c25	= [n16; n19];
c26	= [n19; n1];
c27	= [n19; n11];
c28	= [n5; n18];
c29	= [n18; n21];
c30	= [n21; n5];
c31	= [n21; n15];
c32	= [n18; n9];
c33	= [n9; n20];
c34	= [n20; n13];
c35	= [n16; n5];

% Scale cables by compiling as matrix and scaling matrix
cableMatrix = scalingFactor*[c0; c1; c2; c3; c4; c5; c6; c7; c8; c9; c10;
    c11; c12; c13; c14; c15; c16; c17; c18; c19; c20; c21; c22; c23; c24;
    c25; c26; c27; c28; c29; c30; c31; c32; c33; c34; c35];

% Add cables to plot
for i = 1:2:71
    hold on
    plot3(cableMatrix(i:i+1,1),cableMatrix(i:i+1,2),cableMatrix(i:i+1,3),'b','LineWidth',1.25)
end