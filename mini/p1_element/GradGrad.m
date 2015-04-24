function [A,JF] = GradGrad(P,nu,E)
% Jacobian of the Trasformation
DF(1,1) = P(3,1)-P(1,1);
DF(1,2) = P(5,1)-P(1,1);
DF(2,1) = P(4,1)-P(2,1);
DF(2,2) = P(6,1)-P(2,1);
%
JF = DF(1,1)*DF(2,2)-DF(1,2)*DF(2,1);
% Inverse of Jacobian Trasformation
DFF(1,1) = DF(2,2)/JF;
DFF(1,2) = -DF(2,1)/JF;
DFF(2,1) = -DF(1,2)/JF;
DFF(2,2) = DF(1,1)/JF;
% - Shape functions 1-6    -----------------------------------------------%
epsi(:,1) = DFF*[-1; -1];
epsi(:,2) = DFF*[1; 0];
epsi(:,3) = DFF*[0; 1];
% ------------------------------------------------------------------------%
gradu(:,1) = [epsi(1,1); 0; epsi(2,1)]; % epsiXX, epsiYY, gammaXY
gradu(:,2) = [0; epsi(2,1); epsi(1,1)];
gradu(:,3) = [epsi(1,2); 0; epsi(2,2)];
gradu(:,4) = [0; epsi(2,2); epsi(1,2)];
gradu(:,5) = [epsi(1,3); 0; epsi(2,3)];
gradu(:,6) = [0; epsi(2,3); epsi(1,3)];
% - Legame    ------------------------------------------------------------%
lambda = E*nu/( (1+nu)*(1-2*nu) ); 
G = E/(2*(1+nu));
C = [lambda+2*G, lambda, 0;lambda, lambda+2*G, 0; 0, 0, G];

A = zeros(6,6);
for i =1:6
    for j=1:6
        sigu = C*gradu(:,i);
        A(i,j)=(sigu(1,1)*gradu(1,j)+sigu(2,1)*gradu(2,j)...
            + sigu(3,1)*gradu(3,j))*JF/2;
    end
end
return