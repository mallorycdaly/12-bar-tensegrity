%% Simple example of force-density approach
% Simple, 2D example for experimentation with Schek's force-density
% approach.
%
% See H.-J. Schek's "The Force Density Method for Form Finding and
% Computation of General Networks" for explanation of connectivity matrices
clear;

%% Example from Fig. 2 of paper

% Free nodes connectivity matrix
C = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0;
     1 0 0 -1 0;
     1 -1 0 0 0;
     0 1 -1 0 0;
     0 0 1 -1 0;
     0 0 0 1 -1;
     1 0 0 0 -1;
     0 1 0 0 -1;
     0 0 1 0 -1];
 
% Fixed nodes connectivity matrix
Cf = zeros(12,4);
Cf(1:4,1:4) = -eye(4);
 
% Full connectivity matrix
Cs = [C Cf];

% Nodal positions of fixed nodes defined arbitrarily
xf = [-2 0 2 0]';
yf = [0 -2 0 2]';

% % Coordinate differences
% % Vector form
% u = C*x + Cf*xf;
% v = C*y + Cf*yf;
% % Matrix form
% U = diag(u);
% V = diag(v);

% Select force densities
% q = [ones(4,1); -ones(4,1); ones(4,1)];
q = -ones(12,1);
Q = diag(q);

% Select external forces on free nodes
px = 2*[0.25 0 -0.25 0 0]';
py = px;

% Define D and Df
D = C'*Q*C;
Df = C'*Q*Cf;

% Find free positions
x = inv(D)*(px - Df*xf)
y = inv(D)*(py - Df*yf)

% Plot results
figure
plot(xf,yf,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',7)
hold on
plot(x,y,'o','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',7)
legend('Fixed nodes','Free nodes')

%% Symbolic evaluation

% Symbolic vectors for free positions of nodes
x_sym = sym('x_sym', [1 5], 'real')';
y_sym = sym('y_sym', [1 5], 'real')';

% Symbolic vectors for fixed positions of nodes
xf_sym = sym('x_sym', [1 4], 'real')';
yf_sym = sym('y_sym', [1 4], 'real')';

% Symbolic vector for the force densities (b = bar, c = cable)
q_sym = sym('q', [1 12], 'real')';
Q_sym = diag(q_sym);

% Symbolic coordinate differences
% Vector form
u_sym = C*x_sym + Cf*xf_sym;
v_sym = C*y_sym + Cf*yf_sym;
% Matrix form
U_sym = diag(u_sym);
V_sym = diag(v_sym);

% Equilbrium equations
px_sym = simplify(C'*Q_sym*C*x_sym + C'*Q_sym*Cf*xf_sym);
py_sym = simplify(C'*Q_sym*C*y_sym + C'*Q_sym*Cf*yf_sym);

% Define D and Df for simplicity
D_sym = C'*Q_sym*C;
Df_sym = C'*Q_sym*Cf;