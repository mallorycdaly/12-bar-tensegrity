%% Check Cube Connectivity Matrix
% Created: 2/23/17
% Last Updated: 2/23/17

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
    c11; c12; c13; c14; c15; c16; c17; c18; c19; c20; c21; c22; c23; c24;
    c25; c26; c27; c28; c29; c30; c31; c32; c33; c34; c35] + 1;
% Note: added 1 to change to Matlab indexing

% Number of cables and nodes in 12-bar tensegrity system
numCables = 36;
numNodes = 24;

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