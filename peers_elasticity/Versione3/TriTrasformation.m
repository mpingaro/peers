function [A,DF,DFF,JF] = TriTrasformation(P)
%   P vector 6x1 include the vertex coordinates of the triangle
%   P = [P1x P1y P2x P2y P3x P3y]

DF(1,1) = P(3)-P(1);
DF(2,1) = P(5)-P(1);
DF(3,1) = P(4)-P(2);
DF(4,1) = P(6)-P(2);
%
JF = DF(1,1)*DF(4,1)-DF(2,1)*DF(3,1);
A = JF/2;
%
DFF(1,1) = DF(4,1)/JF;
DFF(2,1) = -DF(2,1)/JF;
DFF(3,1) = -DF(3,1)/JF;
DFF(4,1) = DF(1,1)/JF;

return