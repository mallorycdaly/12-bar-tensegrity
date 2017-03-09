%% Find Equilibrium Position of 12-Bar Tensegrity Octahedron
% Author: Mallory Daly
% Created: 2/22/17
% Last Updated: 2/22/17
clear; close all;

% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks" for explanation of connectivity matrix

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
%     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
Cc = [1  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  0
      1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0;  %  1
      1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0;  %  2
      0  0  1  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  3
      0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0;  %  4
      0  0  0  0  1  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  5
      0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0;  %  6
      0  0  0  0  0  0  1  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  7
      0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0;  %  8
      0  0  0  0  0  1  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  %  9
      0  0  0  0  0  0  0  0  1  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0;  % 10
      0  0  0  0  0  1  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0;  % 11
      0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0;  % 12
      0  0  0  0  0  0  0  0  0  0  1  0 -1  0  0  0  0  0  0  0  0  0  0  0;  % 13
      0  0  0  0  0  0  0  0  0  0  1  0  0  0  0 -1  0  0  0  0  0  0  0  0;  % 14
      0  0  0  0  0  0  0  0  0  0  0  0  1  0 -1  0  0  0  0  0  0  0  0  0;  % 15
      0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0 -1;  % 16
      0  1  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0;  % 17
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0 -1  0  0  0  0  0  0  0;  % 18
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0 -1  0  0  0  0  0;  % 19
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0 -1;  % 20
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0 -1  0  0  0;  % 21
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0 -1  0;  % 22
      0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0 -1  0  0  0;  % 23
      0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -1  0  0;  % 24
      0  0  0  1  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0;  % 25
      0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0 -1  0  0  0  0  0  0;  % 26
      0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0;  % 27
      0  1  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  % 28
      0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0;  % 29
      0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0;  % 30
      0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0 -1  0  0  0  0;  % 31
      0  0  0  0  0  0  0  1  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0;  % 32
      0  0  0  1  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  % 33
      0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0 -1;  % 34
      0  0  0  0  0  0  0  0  0  1  0  0  0 -1  0  0  0  0  0  0  0  0  0  0]; % 35

% Check that there are three entries per node in the cable connectivity
% matrix
nodeEntries = zeros(1,24);
for n = 1:24
    for m = 1:36
        if Cc(m,n) ~= 0
            nodeEntries(n) = nodeEntries(n) + 1;
        end
    end
end
if all(nodeEntries == 3)
    disp('Three entries per node in the cable connectivity matrix!')
end

% Check that there are two entries per cable in the cable connectivity
% matrix
cableEntries = zeros(1,36);
for m = 1:36
    for n = 1:24
        if Cc(m,n) ~= 0
            cableEntries(m) = cableEntries(m) + 1;
        end
    end
end
if all(cableEntries == 2)
    disp('Two entries per cable in the cable connectivity matrix!')
end

% Full Connectivity Matrix
C = [Cb; Cc];

%% Symbolic Evaluation

% Symbolic arrays for x,y,z positions of nodes
x = [sym('x0','real') sym('x', [1 23], 'real')]';
y = [sym('y0','real') sym('y', [1 23], 'real')]';
z = [sym('z0','real') sym('z', [1 23], 'real')]';

% Symbolic arrays for the force densities (b = bar, c = cable)
qb = [sym('qb0','real') sym('qb', [1 11], 'real')]'; % enforce constraint that these are >0
qc = [sym('qc0','real') sym('qc', [1 35], 'real')]'; % enforce constraint that these are <0
q = [qb; qc];

% Symbolic equations
px = C'*diag(C*x)*q;
py = C'*diag(C*y)*q;
pz = C'*diag(C*z)*q;

%% Initial Guess at Node Positions
% From SW model 'C:\Users\Mallory\Google Drive\Tensegrity Research at NASA
% Ames and UC Berkeley\12-Bar Tensegrity\12-Bar Structures\Octahedron\
% octahedron adjustable dims.sldprt'
% See also 'C:\Users\Mallory\Google Drive\Tensegrity Research at NASA Ames
% and UC Berkeley\12-Bar Tensegrity\12-Bar Structures\Node Finding\node
% numbers.ppx' for numbering system
n0  = [-33.37 5.83 5.82];
n1  = [-11.63 -31.81 -5.82];
n2  = [-19.49 25.98 -11.25];
n3  = [19.49 25.98 11.25];
n4  = [-2.53 30.36 -15.91];
n5  = [19.2 -7.28 -27.56];
n6  = [-3.3 20.27 -27.56];
n7  = [-25.03 -17.37 -15.91];
n8  = [22.5 12.99 -22.5];
n9  = [22.5 12.99 22.5];
n10 = [30.74 -12.99 -8.24];
n11 = [-8.24 -12.99 -30.74];
n12 = [25.03 -17.37 15.91];
n13 = [3.30 20.27 27.56];
n14 = [11.63 -31.81 5.82];
n15 = [33.37 5.83 -5.82];
n16 = [0 -25.98 22.5];
n17 = [0 -25.98 -22.5];
n18 = [-19.2 -7.28 27.56];
n19 = [2.53 30.36 15.91];
n20 = [-22.5 12.99 22.5];
n21 = [-22.5 12.99 -22.5];
n22 = [-30.74 -12.99 8.24];
n23 = [8.24 -12.99 30.74];

% Scale nodes by compiling as matrix and scaling matrix
scalingFactor = 0.1;
nMatrix = scalingFactor*[n0; n1; n2; n3; n4; n5; n6; n7; n8; n9; n10; n11; 
    n12; n13; n14; n15; n16; n17; n18; n19; n20; n21; n22; n23];

% Extract initial nodal positions in x,y,z coordinates
xInit = nMatrix(:,1);
yInit = nMatrix(:,2);
zInit = nMatrix(:,3);

% Rod length (45 cm)
rodLength = scalingFactor*45;

%% Find Forces Given Node Positions

% Minimum Tension
minTension = 1;

% Solve with YALMIP
yalmip('clear')
q = sdpvar(48,1);
obj = q'*q;
constr = [C'*diag(C*xInit)*q == 0;
          C'*diag(C*yInit)*q == 0;
          C'*diag(C*zInit)*q == 0;
          q(1:12) >= 0;
          q(13:48) <= -minTension];
options = sdpsettings('solver','quadprog','verbose',1);
sol = optimize(constr,obj,options);

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