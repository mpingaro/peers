function [COORDINATES] = coordcook(NODES,NX,NY,DL1,DL2)
%% THETA
THETA = zeros(NY+1,1);
for I = 0:length(THETA)-1
    THETA(I+1,1) = atand((NODES(2,2)+I*DL1/NY-I*DL2/NY)/NODES(2,1));
end
%% COORDINATES MATRIX
% This matrix contains at rows the nodal points of physical element and at 
% columns the relative coordinates
COORDINATES = zeros((NX+1)*(NY+1),2);
for I = 1:NY+1
    for J = 1:NX+1
        COORDINATES(J+(I-1)*(NX+1),1) = (J-1)*NODES(2,1)/NX;
        COORDINATES(J+(I-1)*(NX+1),2) = COORDINATES(J+(I-1)*(NX+1),1)*tand(THETA(I,1))+...
            (I-1)*NODES(4,2)/NY;
    end
end
end

