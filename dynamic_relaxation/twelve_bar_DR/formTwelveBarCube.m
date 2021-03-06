function [r, cable_pair, rod_pair, L_cable, L_rod] = ...
    formTwelveBarCube(extra_cables)
% This function forms a 12-bar tensegrity cube. Node positions found using
% a SolidWorks model of the structure. 
%
% The inputs are the following:
%   scaling_factor: multiplicative scaling of nodes (0.1 imitates NTRT)
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
scaling_factor = 0.1;
n0  = [-20.79 20.79 8.61]*scaling_factor;
n1  = [-8.61 -20.79 20.79]*scaling_factor;
n2  = [-20.79 20.79 -8.61]*scaling_factor;
n3  = [20.79 8.61 -20.79]*scaling_factor;
n4  = [-8.61 20.79 -20.79]*scaling_factor;
n5  = [-20.79 -20.79 -8.61]*scaling_factor;
n6  = [8.61 20.79 -20.79]*scaling_factor;
n7  = [20.79 8.61 20.79]*scaling_factor;
n8  = [20.79 20.79 -8.61]*scaling_factor;
n9  = [8.61 -20.79 -20.79]*scaling_factor;
n10 = [20.79 20.79 8.61]*scaling_factor;
n11 = [-20.79 8.61 20.79]*scaling_factor;
n12 = [8.61 20.79 20.79]*scaling_factor;
n13 = [20.79 -20.79 8.61]*scaling_factor;
n14 = [-8.61 20.79 20.79]*scaling_factor;
n15 = [-20.79 8.61 -20.79]*scaling_factor;
n16 = [-20.79 -20.79 8.61]*scaling_factor;
n17 = [20.79 -8.61 20.79]*scaling_factor;
n18 = [-8.61 -20.79 -20.79]*scaling_factor;
n19 = [-20.79 -8.61 20.79]*scaling_factor;
n20 = [20.79 -20.79 -8.61]*scaling_factor;
n21 = [-20.79 -8.61 -20.79]*scaling_factor;
n22 = [8.61 -20.79 20.79]*scaling_factor;
n23 = [20.79 -8.61 -20.79]*scaling_factor;
r = [n0; n1; n2; n3; n4; n5; n6; n7; n8; n9; n10; n11; n12; n13; n14; ...
     n15; n16; n17; n18; n19; n20; n21; n22; n23];
 
% % More spherical cube (triangle segments are larger than connecting
% % segments at ~2:1 ratio)
% r = [ ...   
%    -1.8027    1.9885    0.4612;
%    -0.4612   -1.9885    1.8027;
%    -1.9885    1.8027   -0.4612;
%     1.9885    0.4612   -1.8027;
%    -0.4612    1.9885   -1.8027;
%    -1.8027   -1.9885   -0.4612;
%     0.4612    1.8027   -1.9885;
%     1.8027    0.4612    1.9885;
%     1.8027    1.9885   -0.4612;
%     0.4612   -1.9885   -1.8027;
%     1.9885    1.8027    0.4612;
%    -1.9885    0.4612    1.8027;
%     0.4612    1.9885    1.8027;
%     1.8027   -1.9885    0.4612;
%    -0.4612    1.8027    1.9885;
%    -1.8027    0.4612   -1.9885;
%    -1.9885   -1.8027    0.4612;
%     1.9885   -0.4612    1.8027;
%    -0.4612   -1.8027   -1.9885;
%    -1.8027   -0.4612    1.9885;
%     1.9885   -1.8027   -0.4612;
%    -1.9885   -0.4612   -1.8027;
%     0.4612   -1.8027    1.9885;
%     1.8027   -0.4612   -1.9885];

% % Final configuration from octagon-triangle step to land on ground
% % defined by nodes [9 20 23]
% r = [ ...
%    -1.5076    1.5330    0.3082
%    -0.4348   -2.5904    1.7135
%    -2.1672    1.3809   -1.2857
%     2.2737    1.9146   -1.2158
%    -0.7497    1.8484   -2.1906
%    -2.2671   -2.1812   -0.9573
%     0.9573    2.2671   -2.1812
%     2.1906    0.7497    1.8484
%     2.4528    2.2606   -1.2144
%    -0.3749   -1.1299   -0.9686
%     2.5904    1.7135    0.4348
%    -1.5330    0.3082    1.5076
%     1.2857    2.1672    1.3809
%     1.2158   -2.2737    1.9146
%    -0.3082    1.5076    1.5330
%    -1.7135    0.4348   -2.5904
%    -1.8484   -2.1906    0.7497
%     2.1812   -0.9573    2.2671
%    -1.9146   -1.2158   -2.2737
%    -1.3809   -1.2857    2.1672
%     1.1299   -0.9686    0.3749
%    -2.2606   -1.2144   -2.4528
%     1.2144   -2.4528    2.2606
%     0.9686    0.3749   -1.1299];

% Cables by pairs of node indices; add 1 for MATLAB indexing
lattice_cables = [ 0  2;  %  0
                          0 14;  %  1
                          0 11;  %  2
                          2  4;  %  3
                          2 15;  %  4
                          4 15;  %  5
                          4  6;  %  6
                          6  8;  %  7
                          6  3;  %  8
                          8  3;  %  9
                          3 23;  % 10
                          8 10;  % 11
                         10 12;  % 12
                         10  7;  % 13
                          7 12;  % 14
                         14 12;  % 15
                         14 11;  % 16
                         23 20;  % 17
                          9 23;  % 18
                          7 17;  % 19
                         17 13;  % 20
                         13 22;  % 21
                         22 17;  % 22
                         22  1;  % 23
                          1 16;  % 24
                         16 19;  % 25
                         19  1;  % 26
                         19 11;  % 27
                          5 18;  % 28
                         18 21;  % 29
                         21  5;  % 30
                         21 15;  % 31
                         18  9;  % 32
                          9 20;  % 33
                         20 13;  % 34
                         16  5] + 1;  % 35

% Add actuated cables if they exist                     
if isempty(extra_cables) == 1
    cable_pair = lattice_cables;
else
    cable_pair = [lattice_cables; extra_cables];
end
           
% Rods by pairs of node indices; add 1 for MATLAB indexing        
rod_pair   = [ 0  1;
               2  3;
               4  5;
               6  7;
               8  9;
              10 11;
              12 13;
              14 15;
              16 17;
              18 19;
              20 21;
              22 23] + 1;

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

% Ground faces
% Note: Node indices of ground faces should be ordered so that sequential
% indices form an edge. This doesn't matter for a triangle, but matters for
% any larger polygon. Important to the calculation of the distance of the
% projected COG to the edge.
% ground_face = [ 1 16  5 18  9 20 13 22;
%                 0  2 15 21  5 16 19 11;
%                 0  2  4  6  8 10 12 14;
%                10  8  3 23 20 13 17  7;
%                14 12  7 17 22  1 19 11;
%                 4  6  3 23  9 18 21 15] + 1;
% ground_face = [   4  6  3 23  9 18 21 15] + 1;

% Since the triangles won't be the base face, I don't need to include them
% in my analysis            
% ground_face{2} = [ 0 14 11;
%                    2  4 15;
%                    6  8  3;
%                   10 12  7;
%                    1 16 19;
%                    5 21 18;
%                    9 23 20;
%                   13 17 22] + 1;