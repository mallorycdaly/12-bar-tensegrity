function [r, cable_pair, rod_pair, L_cable, L_rod] = ...
    formTwelveBarOctahedron(extra_cables)
% This function forms a 12-bar tensegrity cube. Node positions found using
% a SolidWorks model of the structure. 
%
% The inputs are the following:
%   extra_cables: each row defines the node indices corresponding to a
%       cable added to the base 12-bar structure (with MATLAB indexing)
%
% The outputs are the following:
%   r: 3D array of x,y,z, position of nodes across simulation steps
%   cable_pair: each row defines the node indices corresponding to that
%       cable; the actuated cable pairs are placed at the end
%   rod_pair: each row defines the node indices corresponding to that rod
%   L_cable: vector of cable lengths
%   L_rod: vector of rod length

% Rod length = 45 cm
% scaling_factor = 0.1;
% n0 = [0 0 0]*scaling_factor;
% n1 = [5.4 3.5 18.4]*scaling_factor;
% n2 = [23.2 0 13.8]*scaling_factor;
% n3 = [34.9 3.5 0]*scaling_factor;
% n4 = [23.2 0 -13.8]*scaling_factor;
% n5 = [5.4 3.5 -18.4]*scaling_factor;
% n6 = [-6.2 18 -18.1]*scaling_factor;
% n7 = [-12.2 14 0]*scaling_factor;
% n8 = [9.7 18 27]*scaling_factor;
% n9 = [27.2 14 23.8]*scaling_factor;
% n10 = [41.8 18 -8.5]*scaling_factor;
% n11 = [28.9 14 -22.5]*scaling_factor;
% n12 = [6.8 27 -24.6]*scaling_factor;
% n13 = [-12.4 31 2.5]*scaling_factor;
% n14 = [-4 27 18.2]*scaling_factor;
% n15 = [29.5 31 21.3]*scaling_factor;
% n16 = [39.5 27 4.6]*scaling_factor;
% n17 = [26 31 -24.3]*scaling_factor;
% n18 = [0 41.5 -15.5]*scaling_factor;
% n19 = [0 45 0]*scaling_factor;
% n20 = [3.6 41.5 20.1]*scaling_factor;
% n21 = [23.2 45 13.8]*scaling_factor;
% n22 = [36.5 41.5 0]*scaling_factor;
% n23 = [23.2 45 -13.8]*scaling_factor;
% r = [n0; n1; n2; n3; n4; n5; n6; n7; n8; n9; n10; n11; n12; n13; n14; ...
%      n15; n16; n17; n18; n19; n20; n21; n22; n23];

% Positions found using equilbrium finder with rest length at 95%. All rods
% are the same length (4.4919) whereas original coordinates do not have
% equal rod lengths.
r = [   -0.0340   -0.0465   -0.0385
    0.4888    0.3884    1.6978
    2.2157    0.0640    1.3545
    3.3545    0.4468    0.0622
    2.2954    0.0450   -1.2907
    0.5913    0.3533   -1.7483
   -0.5061    1.7477   -1.8316
   -1.1171    1.3607   -0.1155
    0.9663    1.7678    2.4994
    2.6546    1.4903    2.1578
    3.9612    1.8658   -0.7394
    2.7688    1.4302   -2.0994
    0.6576    2.7199   -2.4475
   -1.1462    2.9616    0.2684
   -0.3156    2.7835    1.7943
    2.9771    3.1201    2.0569
    4.0155    2.7738    0.5031
    2.4838    3.0668   -2.1910
    0.1284    4.0174   -1.3779
   -0.1373    4.4432    0.0602
    0.6272    4.0484    1.8440
    2.5091    4.5444    1.3961
    3.5691    4.0717   -0.0608
    2.3120    4.5354   -1.3542];

% Cables by pairs of node indices; add 1 for MATLAB indexing
lattice_cables = [ 0  1;  %  0
                   0  5;  %  1
                   0  7;  %  2
                   1  2;  %  3
                   1  8;  %  4
                   2  3;  %  5
                   2  9;  %  6
                   3  4;  %  7
                   3 10;  %  8
                   4  5;  %  9
                   4 11;  % 10
                   5  6;  % 11
                   6  7;  % 12
                   6 12;  % 13
                   7 13;  % 14
                   8  9;  % 15
                   8 14;  % 16
                   9 15;  % 17
                  10 11;  % 18
                  10 16;  % 19
                  11 17;  % 20
                  12 17;  % 21
                  12 18;  % 22
                  13 14;  % 23
                  13 19;  % 24
                  14 20;  % 25
                  15 16;  % 26
                  15 21;  % 27
                  16 22;  % 28
                  17 23;  % 29
                  18 19;  % 30
                  18 23;  % 31
                  19 20;  % 32
                  20 21;  % 33
                  21 22;  % 34
                  22 23] + 1;  % 35

% Add actuated cables if they exist                     
if isempty(extra_cables) == 1
    cable_pair = lattice_cables;
else
    cable_pair = [lattice_cables; extra_cables];
end
           
% Rods by pairs of node indices; add 1 for MATLAB indexing        
rod_pair   = [ 0 19;
               2 21;
               4 23;
               6  3;
              13  9;
              18 16;
               1 10;
               7 17;
              14 22;
              20 12;
              15 11;
               8  5] + 1;

% Find number of nodes, cables, and rods           
num_nodes = size(r,1);
num_cables = size(cable_pair,1);
num_actuated_cables = size(extra_cables,1);
num_rods = size(rod_pair,1);

% Original cable and rod lengths, before any deformation
L_cable = zeros(num_cables,1);
for i = 1: num_cables
    L_cable(i) = norm(r(cable_pair(i,1),:) - ...
        r(cable_pair(i,2),:));
end
L_rod = norm(r(rod_pair(1,1),:) - r(rod_pair(1,2),:));