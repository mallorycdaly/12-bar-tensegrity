%% Find Equilibrium Position of 12-Bar Tensegrity Octahedron
% Author: Mallory Daly
% Created: 2/22/17
% Last Updated: 2/22/17
clear; close all;

% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks" for explanation of connectivity matrix

%% System Parameters

% Number of cables and nodes in 12-bar tensegrity system
numCables = 36;
numNodes = 24;

%% Rod Mass from Geometry

% Geometry
OD = 0.375*0.0254; % outer diameter (m)
ID = (0.375-0.035)*0.0254; % inner diameter (m)
L = 0.45; % rod length (m)
V = pi*((OD/2)^2-(ID/2)^2)*L;

% Density of aluminum 6061
p = 2720; % density (kg/m^3)

% Mass per rod
m = p*V; % mass (kg)

%% Connectivity Matrices
% Columns are nodes; rows are bars/cables

% Connectivity matrix for the bars
%     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
Cb = [1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  0
      0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  1
      0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  2
      0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  3
      0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  4
      0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0;  %  5
      0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0;  %  6
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0;  %  7
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0;  %  8
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0;  %  9
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0;  % 10
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1]; % 11

% Connectivity matrix for the cables

    % Cable node connections
    c0	= [0 2];
    c1	= [0 14];
    c2	= [0 11];
    c3	= [2 4];
    c4	= [2 15];
    c5	= [4 15];
    c6	= [4 6];
    c7	= [6 8];
    c8	= [6 3];
    c9	= [8 3];
    c10	= [3 23];
    c11	= [8 10];
    c12	= [10 12];
    c13	= [10 7];
    c14	= [7 12];
    c15	= [14 12];
    c16	= [14 11];
    c17	= [23 20];
    c18	= [9 23];
    c19	= [7 17];
    c20	= [17 13];
    c21	= [13 22];
    c22	= [22 17];
    c23	= [22 1];
    c24	= [1 16];
    c25	= [16 19];
    c26	= [19 1];
    c27	= [19 11];
    c28	= [5 18];
    c29	= [18 21];
    c30	= [21 5];
    c31	= [21 15];
    c32	= [18 9];
    c33	= [9 20];
    c34	= [20 13];
    c35	= [16 5];

    % Matrix of cable node connections
    C = [c0; c1; c2; c3; c4; c5; c6; c7; c8; c9; c10;
        c11; c12; c13; c14; c15; c16; c17; c18; c19; c20;
        c21; c22; c23; c24; c25; c26; c27; c28; c29; c30;
        c31; c32; c33; c34; c35] + 1;
    % note: added 1 to change to Matlab indexing

    % Initialize cable connectivity matrix
    Cc = zeros(numCables,numNodes);

    % Fill in 1 and -1 per cable using cable node connections
    for i = 1:numCables
        if C(i,1) < C(i,2)
            Cc(i,C(i,1)) = 1;
            Cc(i,C(i,2)) = -1;
        else
            Cc(i,C(i,2)) = 1;
            Cc(i,C(i,1)) = -1;
        end
    end

% Full Connectivity Matrix
C = [Cb; Cc];

%% Symbolic Evaluation

% Symbolic arrays for x,y,z positions of nodes
x = [sym('x0','real') sym('x', [1 23], 'real')]';
y = [sym('y0','real') sym('y', [1 23], 'real')]';
z = [sym('z0','real') sym('z', [1 23], 'real')]';

% Symbolic arrays for the force densities (b = bar, c = cable)
% qb = [sym('qb0','real') sym('qb', [1 11], 'real')]'; % enforce constraint that these are >0
% qc = [sym('qc0','real') sym('qc', [1 35], 'real')]'; % enforce constraint that these are <0
% q = [qb; qc];
% q = [sym('q0','real') sym('q', [1 47], 'real')]'; % enforce constraint that these are >0

% Symbolic equations
% px = C'*diag(C*x)*q;
% py = C'*diag(C*y)*q;
% pz = C'*diag(C*z)*q;

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

% Scale nodes by compiling as matrix and scaling matrix
scalingFactor = 0.1;
nMatrix = scalingFactor*[n0; n1; n2; n3; n4; n5; n6; n7; n8; n9; n10; n11; 
    n12; n13; n14; n15; n16; n17; n18; n19; n20; n21; n22; n23];

% Extract initial nodal positions in x,y,z coordinates
xInit = nMatrix(:,1);
yInit = nMatrix(:,2);
zInit = nMatrix(:,3);

% Plot nodes
% figure
% scatter3(xInit,yInit,zInit,'filled')
% title('12-Bar Tensegrity Cube')
% axis equal

% Symbolic equations
% px = C'*diag(C*xInit)*q;
% py = C'*diag(C*yInit)*q;
% pz = C'*diag(C*zInit)*q;

%% Find Cable Tensions and Rod Compressive Forces Given Node Positions

% Minimum tension
minTension = 0;

% Weight per rod
w = 98.1*m; % (kg*dm/s^2) using dm to be consistent with NTRT

% Total system weight
W = 12*w;

% External forces
py = zeros(24,1);
py([1,16,5,18,9,20,13,22]+1) = W/8; % assume only vertical reaction force, distributed evenly among 8 nodes in contact with ground
py = py -w/2; % gravity acts in the negative y direction; mass is distributed equally per node

% Solve with YALMIP
yalmip('clear')
q = sdpvar(48,1);
obj = q'*q;
constr = [C'*diag(C*xInit)*q == 0;
          C'*diag(C*yInit)*q == 0;
%           C'*diag(C*yInit)*q + py == 0; 
          C'*diag(C*zInit)*q == 0;
%           q(1:12) >= 0; % bars can only push
          q(13:48) <= -minTension]; % cables can only pull
options = sdpsettings('solver','quadprog','verbose',1);
% options = sdpsettings('verbose',1);
sol = optimize(constr,obj,options);

% Show solutions
barForceDensities = value(q(1:12))'
cableForceDensities = value(q(13:48))'

%% Find Node Positions Given Forces
% 
% minTension = 10;
% 
% yalmip('clear')
% x = sdpvar(24,1);
% y = sdpvar(24,1);
% z = sdpvar(24,1);
% q = sdpvar(48,1)
% % assign(x,xInit);
% % assign(y,yInit);
% % assign(z,zInit);
% obj = sum((x-xInit).^2+(y-yInit).^2+(z-zInit).^2);
% constr = [C'*diag(C*xInit)*q == 0;
%           C'*diag(C*yInit)*q == 0;
%           C'*diag(C*zInit)*q == 0;
%           q(1:12) >= 0;
%           q(13:48) <= -minTension];
% % options = sdpsettings('verbose',1,'usex0',1);
% options = sdpsettings('verbose',1);
% sol = optimize(constr,obj,options);