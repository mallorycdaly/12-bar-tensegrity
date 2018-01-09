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

for i = 1:length(r)
    fprintf(['n' num2str(i-1) ': [' num2str(r(i,1)) ',' num2str(r(i,2)) ',' num2str(r(i,3)) ']\n'])
end
% 
% scaling_factor = 0.1;
% n0  = [-20.79 20.79 8.61]*scaling_factor;
% n1  = [-8.61 -20.79 20.79]*scaling_factor;
% n2  = [-20.79 20.79 -8.61]*scaling_factor;
% n3  = [20.79 8.61 -20.79]*scaling_factor;
% n4  = [-8.61 20.79 -20.79]*scaling_factor;
% n5  = [-20.79 -20.79 -8.61]*scaling_factor;
% n6  = [8.61 20.79 -20.79]*scaling_factor;
% n7  = [20.79 8.61 20.79]*scaling_factor;
% n8  = [20.79 20.79 -8.61]*scaling_factor;
% n9  = [8.61 -20.79 -20.79]*scaling_factor;
% n10 = [20.79 20.79 8.61]*scaling_factor;
% n11 = [-20.79 8.61 20.79]*scaling_factor;
% n12 = [8.61 20.79 20.79]*scaling_factor;
% n13 = [20.79 -20.79 8.61]*scaling_factor;
% n14 = [-8.61 20.79 20.79]*scaling_factor;
% n15 = [-20.79 8.61 -20.79]*scaling_factor;
% n16 = [-20.79 -20.79 8.61]*scaling_factor;
% n17 = [20.79 -8.61 20.79]*scaling_factor;
% n18 = [-8.61 -20.79 -20.79]*scaling_factor;
% n19 = [-20.79 -8.61 20.79]*scaling_factor;
% n20 = [20.79 -20.79 -8.61]*scaling_factor;
% n21 = [-20.79 -8.61 -20.79]*scaling_factor;
% n22 = [8.61 -20.79 20.79]*scaling_factor;
% n23 = [20.79 -8.61 -20.79]*scaling_factor;
% r = [n0; n1; n2; n3; n4; n5; n6; n7; n8; n9; n10; n11; n12; n13; n14; ...
%      n15; n16; n17; n18; n19; n20; n21; n22; n23]