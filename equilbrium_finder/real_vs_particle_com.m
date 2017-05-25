%% Real vs. particle model to find center of mass
% Show that total center of mass (com) can be found using com of each rod,
% or equally weighted particles placed on ends of rod
clear; close all

%% Original structure
% 3 bar 2D tetrahedra

% Center of mass (com) of each rod
rodMass = 1;
com1 = [0.5; 0.5];
com2 = [1.5; 0.5];
com3 = [1; 1.5];
com_real = 1/(3*rodMass)*(com1*rodMass + com2*rodMass + com3*rodMass)

% Particle approximation of rods (one particle at each end of rod)
particleMass = rodMass/2;
part1 = [0; 0];
part2 = [1; 1];
part3 = [1; 1];
part4 = [1; 1];
part5 = [1; 2];
part6 = [2; 0];
com_particle = 1/(6*rodMass/2)*particleMass*(part1 + part2 + part3 + ...
    part4 + part5 + part6)

%% Rotated structure

% Rotation angle in degrees
rotAngle = 45;
rotM = [cosd(rotAngle) -sind(rotAngle); sind(rotAngle) cosd(rotAngle)];

% Center of mass (com) of each rod
rodMass = 1;
com1 = rotM*[0.5; 0.5];
com2 = rotM*[1.5; 0.5];
com3 = rotM*[1; 1.5];
com_real = 1/(3*rodMass)*(com1*rodMass + com2*rodMass + com3*rodMass)

% Particle approximation of rods (one particle at each end of rod)
particleMass = rodMass/2;
part1 = rotM*[0; 0];
part2 = rotM*[1; 1];
part3 = rotM*[1; 1];
part4 = rotM*[1; 1];
part5 = rotM*[1; 2];
part6 = rotM*[2; 0];
com_particle = 1/(6*rodMass/2)*particleMass*(part1 + part2 + part3 + ...
    part4 + part5 + part6)