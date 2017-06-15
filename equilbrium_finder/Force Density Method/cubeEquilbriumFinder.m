%% Find Equilibrium Position of 12-Bar Tensegrity Octahedron
% Author: Mallory Daly
% Created: 2/22/17
% Last Updated: 5/24/17
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks" for explanation of connectivity matrices

clear; close all;

%% System Parameters

% Number of cables and nodes in 12-bar tensegrity system
numRods = 12;
numCables = 36;
numNodes = 24;

%% Rod Mass from Geometry

% Geometry
rodOuterDiameter = 0.375*0.0254; % (m)
rodInnerDiameter = (0.375-0.035)*0.0254; % (m)
rodLength = 0.45; % (m)
rodVolume = pi*((rodOuterDiameter/2)^2-(rodInnerDiameter/2)^2)*rodLength; % (m^3)

% Density of aluminium 6061
densityAluminium = 2720; % (kg/m^3)

% Mass per rod
rodMass = densityAluminium*rodVolume; % (kg)

%% Connectivity Matrices

% Connectivity matrix for the bars
% Nodes 0 to 1 form bar 0, nodes 2 to 3 form bar 1, and so on
%                                    NODES
%          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
conBars = [1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  0
           0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  1
           0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  2
           0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  3
           0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  4  B
           0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0;  %  5  A
           0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0;  %  6  R
           0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0;  %  7  S
           0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0;  %  8
           0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0;  %  9
           0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0;  % 10
           0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1]; % 11
     
% Connectivity matrix for the cables

    % Node numbers defining each cable
    cable0	= [0 2];
    cable1	= [0 14];
    cable2	= [0 11];
    cable3	= [2 4];
    cable4	= [2 15];
    cable5	= [4 15];
    cable6	= [4 6];
    cable7	= [6 8];
    cable8	= [6 3];
    cable9	= [8 3];
    cable10	= [3 23];
    cable11	= [8 10];
    cable12	= [10 12];
    cable13	= [10 7];
    cable14	= [7 12];
    cable15	= [14 12];
    cable16	= [14 11];
    cable17	= [23 20];
    cable18	= [9 23];
    cable19	= [7 17];
    cable20	= [17 13];
    cable21	= [13 22];
    cable22	= [22 17];
    cable23	= [22 1];
    cable24	= [1 16];
    cable25	= [16 19];
    cable26	= [19 1];
    cable27	= [19 11];
    cable28	= [5 18];
    cable29	= [18 21];
    cable30	= [21 5];
    cable31	= [21 15];
    cable32	= [18 9];
    cable33	= [9 20];
    cable34	= [20 13];
    cable35	= [16 5];

    % Form matrix of cable node numbers and add 1 for Matlab indexing
    cables = [cable0; cable1; cable2; cable3; cable4; cable5; cable6;
        cable7; cable8; cable9; cable10; cable11; cable12; cable13;
        cable14; cable15; cable16; cable17; cable18; cable19; cable20;
        cable21; cable22; cable23; cable24; cable25; cable26; cable27;
        cable28; cable29; cable30; cable31; cable32; cable33; cable34;
        cable35] + 1;

    % Initialize cable connectivity matrix
    conCables = zeros(numCables,numNodes);

    % Fill in 1 and -1 per cable using cable node connections; lower number
    % gets 1 and higher number gets -1 in column position based on node
    % number
    for i = 1:numCables
        if cables(i,1) < cables(i,2)
            conCables(i,cables(i,1)) = 1;
            conCables(i,cables(i,2)) = -1;
        else
            conCables(i,cables(i,2)) = 1;
            conCables(i,cables(i,1)) = -1;
        end
    end

% Full connectivity matrix includes bars and cables
conFull = [conBars; conCables];

%% Symbolic Evaluation

% Symbolic arrays for x,y,z positions of nodes
x = [sym('x0','real') sym('x', [1 23], 'real')]';
y = [sym('y0','real') sym('y', [1 23], 'real')]';
z = [sym('z0','real') sym('z', [1 23], 'real')]';

% % Symbolic arrays for the force densities (b = bar, c = cable)
qb = [sym('qb0','real') sym('qb', [1 11], 'real')]'; % enforce constraint that these are >0
qc = [sym('qc0','real') sym('qc', [1 35], 'real')]'; % enforce constraint that these are <0
q = [qb; qc];
q = [sym('q0','real') sym('q', [1 47], 'real')]'; % enforce constraint that these are >0

% Symbolic equations
px = conFull'*diag(conFull*x)*q;
py = conFull'*diag(conFull*y)*q;
pz = conFull'*diag(conFull*z)*q;

%% Initial Guess at Node Positions
% From SW model 'C:\Users\Mallory\Google Drive\Tensegrity Research at NASA
% Ames and UC Berkeley\12-Bar Tensegrity\12-Bar Structures\Cube\cube
% adjustable dims.sldprt'
% See also 'C:\Users\Mallory\Google Drive\Tensegrity Research at NASA Ames
% and UC Berkeley\12-Bar Tensegrity\12-Bar Structures\Node Finding\node
% numbers.ppx' for numbering system
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

% Scale nodes to duplicate NTRT conditions
scalingFactor = 0.1;
nodes = scalingFactor*[n0; n1; n2; n3; n4; n5; n6; n7; n8; n9; n10; n11; 
    n12; n13; n14; n15; n16; n17; n18; n19; n20; n21; n22; n23];

% Extract initial nodal positions in x,y,z coordinates
xExpected = nodes(:,1);
yExpected = nodes(:,2);
zExpected = nodes(:,3);

% % Plot nodes
% figure
% scatter3(xInit,yInit,zInit,'filled')
% title('12-Bar Tensegrity Cube')
% axis equal
% 
% % Symbolic equations
% px = conFull'*diag(conFull*xInit)*q;
% py = conFull'*diag(conFull*yInit)*q;
% pz = conFull'*diag(conFull*zInit)*q;

%% Find force densities

% Weight per rod
weightPerRod = 98.1*rodMass; % (kg*dm/s^2) using dm to be consistent with NTRT

% Total system weight
weightTotal = 12*weightPerRod;

% External forces
extForce_x = zeros(24,1);
extForce_y = zeros(24,1);
extForce_z = zeros(24,1);
% Only in the y direction from normal forces and gravity
% extForce_y([1,16,5,18,9,20,13,22]+1) = weightTotal/8; % assume only vertical reaction force, distributed evenly among 8 nodes in contact with ground
% extForce_y = extForce_y - weightPerRod/2; % gravity acts in the negative y direction; mass is distributed equally per node

% Minimum tension in the cables (N)
minTension = 10;

% Set up yalmip
yalmip('clear')
x = sdpvar(24,1);
y = sdpvar(24,1);
z = sdpvar(24,1);
q = sdpvar(48,1)

% Minimize difference between calculated and expected node positions
obj = (x-xExpected)'*(x-xExpected) + (y-yExpected)'*(y-yExpected) + ...
    (z-zExpected)'*(z-zExpected);
constr = [(conFull'*diag(q)*conFull)*x == extForce_x;
          (conFull'*diag(q)*conFull)*y == extForce_y;
          (conFull'*diag(q)*conFull)*z == extForce_z];
%            q(1:12) >= 0;
%            q(13:48) < 0];
% options = sdpsettings('verbose',1,'usex0',1);
options = sdpsettings('verbose',1);
sol = optimize(constr,obj,options);

% Extract solutions
xFound = value(x);
yFound = value(y);
zFound = value(z);
qFound = value(q);

%% Plot expected vs. found
figure
scatter3(xExpected,yExpected,zExpected,'filled','blue')
hold on
scatter3(xFound,yFound,zFound,'filled','red')
legend('Expected','Found','Location','Best')
xlabel('x')
ylabel('y')
zlabel('z')
title('12-Bar Tensegrity Cube')
axis equal

%% Find Cable Tensions and Rod Compressive Forces Given Node Positions
% DOESN'T WORK! CAN'T IMPLEMENT MIN TENSION
% 
% % Minimum tension
% minTension = 1;
% 

% 
% % External forces in the y direction from normal forces and gravity
% extForce_y = zeros(24,1);
% extForce_y([1,16,5,18,9,20,13,22]+1) = weightTotal/8; % assume only vertical reaction force, distributed evenly among 8 nodes in contact with ground
% extForce_y = extForce_y - weightPerRod/2; % gravity acts in the negative y direction; mass is distributed equally per node
% %
% % Solve with YALMIP
% yalmip('clear')
% q = sdpvar(48,1); % force density vector of cables and bars
% obj = q'*q;
% constr = [conFull'*diag(conFull*xInit)*q == 0;
%           conFull'*diag(conFull*yInit)*q == 0;
% %           C'*diag(C*yInit)*q + py == 0; 
%           conFull'*diag(conFull*zInit)*q == 0;
% %           q(1:12) >= 0; % bars can only push
%           q(13:48) <= -minTension]; % cables can only pull
% options = sdpsettings('solver','quadprog','verbose',1);
% % options = sdpsettings('verbose',1);
% sol = optimize(constr,obj,options);
% 
% % Show solutions
% barForceDensities = value(q(1:12))
% cableForceDensities = value(q(13:48))
