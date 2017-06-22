function addCoordinateSystemToPlot(r, rod_pair, num_rods)
% If no arguments are input, this function adds the world Cartesian
% coordinate system to a plot. If arguments are input, this function adds
% the world Cartesian coordinate system as well as the local Cartesian
% coordinate system on each rod based on the translation:
%   r = x*E_1 + y*E_2 + z*E_3
% 
% The local coordinate system is based on the Fig. 4 of Kyunam Kim's
% "Robust Learning of Tensegrity Robot Control for Locomotion through
% Form-Finding."
% 
% The inputs are the following:
%   r: matrix of x,y,z, position of nodes
%   rod_pair: each row defines the node indices corresponding to that rod
%   num_rods: number of rods

% World basis, Cartesian
I = 0.2*eye(3);
str = 'xyz';
for i = 1:3
    hold on
    plot3([0 I(i,1)],[0 I(i,2)],[0 I(i,3)],'r','LineWidth',1.5)
    text(I(i,1),I(i,2),I(i,3),[' ' str(i) ' '],'Color','red')
end

if nargin == 3
    
    % Reorder nodes by rod pairs and so that higher (z) point is second
    % Odd nodes: lower value of z
    % Even nodes: higher value of z
    r_rod = zeros(size(r));
    for i = 1:num_rods
        if (r(rod_pair(i,1),3) > r(rod_pair(i,2),3))
            r_rod(2*i-1:2*i,:) = [r(rod_pair(i,2),:)
                r(rod_pair(i,1),:)];
        else
            r_rod(2*i-1:2*i,:) = [r(rod_pair(i,1),:)
                r(rod_pair(i,2),:)];
        end
    end

    % Plot local coordinate systems
    for i = 1:3
        r = r_rod(2*i-1,:);
        for j = 1:3
            hold on
            plot3([0 I(j,1)]+r(1),[0 I(j,2)]+r(2),[0 I(j,3)]+r(3),'r','LineWidth',1)
            text(I(j,1)+r(1),I(j,2)+r(2),I(j,3)+r(3),[' ' str(j) ' '],'Color','red')
        end
        
    end
    
end