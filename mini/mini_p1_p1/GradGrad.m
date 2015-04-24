function [A,B,JF] = GradGrad(P)
%   P vector 6x1 include the vertex coordinates of the triangle
%   P = [P1x P1y P2x P2y P3x P3y]
DF(1,1) = P(3)-P(1);
DF(1,2) = P(5)-P(1);
DF(2,1) = P(4)-P(2);
DF(2,2) = P(6)-P(2);
%
JF = DF(1,1)*DF(2,2)-DF(1,2)*DF(2,1);
%
DFF(1,1) = DF(2,2)/JF;
DFF(1,2) = -DF(2,1)/JF;
DFF(2,1) = -DF(1,2)/JF;
DFF(2,2) = DF(1,1)/JF;

% QUADRATURA DI GAUSS (grado di precisione 5)
a= 9/80;
b=(155+sqrt(15))/2400;
c=(155-sqrt(15))/2400;
w(1,1)=a; w(1,2)=b; w(1,3)=b; w(1,4)=b;
w(1,5)=c; w(1,6)=c; w(1,7)=c;
q1=1/3; q2=(6+sqrt(15))/21; q3=(9-2*sqrt(15))/21;
q4=(6-sqrt(15))/21; q5=(9+2*sqrt(15))/21;
xnod=[q1 q2 q3 q2 q4 q5 q4]; ynod=[q1 q2 q2 q3 q4 q4 q5];

% Shape functions 1-6
epsi(:,1) = DFF*[-1; -1];
epsi(:,2) = DFF*[1; 0];
epsi(:,3) = DFF*[0; 1];

gradu(:,1) = [epsi(1,1); 0; epsi(2,1)/2];
gradu(:,2) = [0; epsi(2,1); epsi(1,1)/2];
gradu(:,3) = [epsi(1,2); 0; epsi(2,2)/2];
gradu(:,4) = [0; epsi(2,2); epsi(1,2)/2];
gradu(:,5) = [epsi(1,3); 0; epsi(2,3)/2];
gradu(:,6) = [0; epsi(2,3); epsi(1,3)/2];
%
A = zeros(8,8);
B = zeros(8,3);
for k=1:7
    x = xnod(1,k); y = ynod(1,k);
    % Shape Functions 7 & 8
    epsi(:,4) = DFF*[27*(y-2*x*y-y*y); 27*(x-2*x*y-x*x)];
    gradu(:,7) = [epsi(1,4); 0; epsi(2,4)/2];
    gradu(:,8) = [0; epsi(2,4); epsi(1,4)/2];
    
    p(1,1) = 1-x-y; p(2,1) = x; p(3,1) = y;

    for i =1:8
        for j=1:8
            A(i,j) = A(i,j) + w(1,k)*( gradu(1,i)*gradu(1,j) +...
                gradu(2,i)*gradu(2,j) + 2*gradu(3,i)*gradu(3,j))*JF;
        end
        for ii =1:3
            divu = gradu(1,i)+gradu(2,i);
            B(i,ii) = B(i,ii) + w(1,k)*p(ii,1)*divu*JF;
        end
    end
end
return