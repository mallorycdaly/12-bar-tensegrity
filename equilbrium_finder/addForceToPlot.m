function addForceToPlot(r, F, S)
% This function adds one or more force vectors to a plot.
%
% The inputs are the following:
%   r: matrix of nodal x,y,z positions
%   F: matrix of force vectors
%   S: character string made from elements of plotting columns (e.g. '-k')

if nargin == 2

    for i = 1:size(r,1)
        x = [0; F(i,1)]+r(i,1);
        y = [0; F(i,2)]+r(i,2);
        z = [0; F(i,3)]+r(i,3);
        hold on
        plot3(x,y,z,'r')
    end
    
elseif nargin == 3
    
    for i = 1:size(r,1)
        x = [0; F(i,1)]+r(i,1);
        y = [0; F(i,2)]+r(i,2);
        z = [0; F(i,3)]+r(i,3);
        hold on
        plot3(x,y,z,S)
    end
    
end