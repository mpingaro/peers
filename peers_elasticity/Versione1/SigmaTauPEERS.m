function [A,B,C,AT] = SigmaTauPEERS(P,mu,lambda,s)
% -- Computing the bilinear form a(Sigma,Tau) e b(u,div Tau)
%   P vector 6x1 include the vertex coordinates of the triangle
%   P = [P1x P1y P2x P2y P3x P3y]
%   mu, lambda = costanti di Lam�
%   s = segno delle funzioni di forma (inversione flussi)
F(1,1) = P(3,1)-P(1,1);
F(1,2) = P(5,1)-P(1,1);
F(2,1) = P(4,1)-P(2,1);
F(2,2) = P(6,1)-P(2,1);
%
FF = inv(F);
dF = F(1,1)*F(2,2)-F(1,2)*F(2,1);
AT = dF/2;
% QUADRATURA DI GAUSS (grado di precisione 5)
nquad = 7;
a= 9/80;
b=(155+sqrt(15))/2400;
c=(155-sqrt(15))/2400;
w(1)=a; w(2)=b; w(3)=b; w(4)=b;
w(5)=c; w(6)=c; w(7)=c;
q1=1/3; q2=(6+sqrt(15))/21; q3=(9-2*sqrt(15))/21;
q4=(6-sqrt(15))/21; q5=(9+2*sqrt(15))/21;
xnod=[q1 q2 q3 q2 q4 q5 q4]; ynod=[q1 q2 q2 q3 q4 q4 q5];
%
% Elastic constant
t1 = 1/(2*mu); t2 = lambda/(4*mu*(mu+lambda));
%
A = zeros(8,8);    % Storo solo triangolare superiore
C = zeros(8,3);    % Storo completamente la matrice
for k=1:nquad
    x(1) = xnod(k); x(2) = ynod(k);
    % Moltiplicatori recupero simmetria
    eta(1)= 1-x(1)-x(2);
    eta(2)= x(1);
    eta(3)= x(2);
    % RT0
    psi(:,1) = F*[x(1); -1+x(2)]/dF;
    psi(:,2) = F*[(2/sqrt(2))*x(1); (2/sqrt(2))*x(2)]/dF;
    psi(:,3) = F*[-1+x(1); x(2)]/dF;
    % BOLLA
    dBx = (x(2)-x(2)*x(2)-2*x(1)*x(2))*27; 
    dBy = (x(1)-x(1)*x(1)-2*x(1)*x(2))*27;
    dBT = FF'*[dBx;dBy];
    psi(:,4) = [dBT(2); -dBT(1)]*dF; 
    % MATRIX A a( C*sigma,tau )
    A(1,1) = A(1,1)+w(k)*(t1*(psi(1,1)*psi(1,1)...
        +psi(2,1)*psi(2,1))-t2*psi(1,1)*psi(1,1))*dF;        
    A(1,2) = A(1,2)+w(k)*(-t2*psi(1,1)*psi(2,1))*dF;
    A(1,3) = A(1,3)+w(k)*(t1*(psi(1,1)*psi(1,2)...
        +psi(2,1)*psi(2,2))-t2*psi(1,1)*psi(1,2))*dF;
    A(1,4) = A(1,4)+w(k)*(-t2*psi(1,1)*psi(2,2))*dF;
    A(1,5) = A(1,5)+w(k)*(t1*(psi(1,1)*psi(1,3)...
        +psi(2,1)*psi(2,3))-t2*psi(1,1)*psi(1,3))*dF;
    A(1,6) = A(1,6)+w(k)*(-t2*psi(1,1)*psi(2,3))*dF;
    %
    A(2,2) = A(2,2)+w(k)*(t1*(psi(1,1)*psi(1,1)...
        +psi(2,1)*psi(2,1))-t2*psi(2,1)*psi(2,1))*dF;
    A(2,3) = A(2,3)+w(k)*(-t2*psi(2,1)*psi(1,2))*dF;
    A(2,4) = A(2,4)+w(k)*(t1*(psi(1,1)*psi(1,2)...
        +psi(2,1)*psi(2,2))-t2*psi(2,1)*psi(2,2))*dF;
    A(2,5) = A(2,5)+w(k)*(-t2*psi(2,1)*psi(1,3))*dF;
    A(2,6) = A(2,6)+w(k)*(t1*(psi(1,1)*psi(1,3)...
        +psi(2,1)*psi(2,3))-t2*psi(2,1)*psi(2,3))*dF;
    %
    A(3,3) = A(3,3)+w(k)*(t1*(psi(1,2)*psi(1,2)...
        +psi(2,2)*psi(2,2))-t2*psi(1,2)*psi(1,2))*dF;
    A(3,4) = A(3,4)+w(k)*(-t2*psi(1,2)*psi(2,2))*dF;
    A(3,5) = A(3,5)+w(k)*(t1*(psi(1,2)*psi(1,3)...
        +psi(2,2)*psi(2,3))-t2*psi(1,2)*psi(1,3))*dF;
    A(3,6) = A(3,6)+w(k)*(-t2*psi(1,2)*psi(2,3))*dF;
    %
    A(4,4) = A(4,4)+w(k)*(t1*(psi(1,2)*psi(1,2)...
        +psi(2,2)*psi(2,2))-t2*psi(2,2)*psi(2,2))*dF;
    A(4,5) = A(4,5)+w(k)*(-t2*psi(2,2)*psi(1,3))*dF;
    A(4,6) = A(4,6)+w(k)*(t1*(psi(1,2)*psi(1,3)...
        +psi(2,2)*psi(2,3))-t2*psi(2,2)*psi(2,3))*dF;
    %
    A(5,5) = A(5,5)+w(k)*(t1*(psi(1,3)*psi(1,3)...
        +psi(2,3)*psi(2,3))-t2*psi(1,3)*psi(1,3))*dF;
    A(5,6) = A(5,6)+w(k)*(-t2*psi(1,3)*psi(2,3))*dF;
    
    A(6,6) = A(6,6)+w(k)*(t1*(psi(1,3)*psi(1,3)...
        +psi(2,3)*psi(2,3))-t2*psi(2,3)*psi(2,3))*dF;
    %
    A(7,7) = A(7,7)+w(k)*(t1*(psi(1,4)*psi(1,4)...
        +psi(2,4)*psi(2,4))-t2*psi(1,4)*psi(1,4))*dF;
    A(7,8) = A(7,8)+w(k)*(-t2*psi(1,4)*psi(2,4))*dF;
    A(8,8) = A(8,8)+w(k)*(t1*(psi(1,4)*psi(1,4)...
        +psi(2,4)*psi(2,4))-t2*psi(2,4)*psi(2,4))*dF;
    %
    A(1,7) = A(1,7)+w(k)*s*(t1*(psi(1,1)*psi(1,4)... 
        +psi(2,1)*psi(2,4))-t2*psi(1,1)*psi(1,4))*dF; 
    A(1,8) = A(1,8)+w(k)*s*(-t2*psi(1,1)*psi(2,4))*dF;
    A(2,7) = A(2,7)+w(k)*s*(-t2*psi(2,1)*psi(1,4))*dF;
    A(2,8) = A(2,8)+w(k)*s*(t1*(psi(1,1)*psi(1,4)...
        +psi(2,1)*psi(2,4))-t2*psi(2,1)*psi(2,4))*dF;
    A(3,7) = A(3,7)+w(k)*s*(t1*(psi(1,2)*psi(1,4)...
        +psi(2,2)*psi(2,4))-t2*psi(1,2)*psi(1,4))*dF;
    A(3,8) = A(3,8)+w(k)*s*(-t2*psi(1,2)*psi(2,4))*dF; 
    A(4,7) = A(4,7)+w(k)*s*(-t2*psi(2,2)*psi(1,4))*dF;
    A(4,8) = A(4,8)+w(k)*s*(t1*(psi(1,2)*psi(1,4)...
        +psi(2,2)*psi(2,4))-t2*psi(2,2)*psi(2,4))*dF;
    A(5,7) = A(5,7)+w(k)*s*(t1*(psi(1,3)*psi(1,4)...
        +psi(2,3)*psi(2,4))-t2*psi(1,3)*psi(1,4))*dF;
    A(5,8) = A(5,8)+w(k)*s*(-t2*psi(1,3)*psi(2,4))*dF;
    A(6,7) = A(6,7)+w(k)*s*(-t2*psi(2,3)*psi(1,4))*dF;
    A(6,8) = A(6,8)+w(k)*s*(t1*(psi(1,3)*psi(1,4)...
        +psi(2,3)*psi(2,4))-t2*psi(2,3)*psi(2,4))*dF;
    % MATRIX C c( eta,as(tau) )
    C(1,1) = C(1,1)+w(k)*s*(-psi(2,1)*eta(1))*dF;
    C(2,1) = C(2,1)+w(k)*s*(psi(1,1)*eta(1))*dF;
    C(3,1) = C(3,1)+w(k)*s*(-psi(2,2)*eta(1))*dF;
    C(4,1) = C(4,1)+w(k)*s*(psi(1,2)*eta(1))*dF;
    C(5,1) = C(5,1)+w(k)*s*(-psi(2,3)*eta(1))*dF;
    C(6,1) = C(6,1)+w(k)*s*(psi(1,3)*eta(1))*dF;
    C(7,1) = C(7,1)+w(k)*(-psi(2,4)*eta(1))*dF;
    C(8,1) = C(8,1)+w(k)*(psi(1,4)*eta(1))*dF;
    %  
    C(1,2) = C(1,2)+w(k)*s*(-psi(2,1)*eta(2))*dF;
    C(2,2) = C(2,2)+w(k)*s*(psi(1,1)*eta(2))*dF;
    C(3,2) = C(3,2)+w(k)*s*(-psi(2,2)*eta(2))*dF;
    C(4,2) = C(4,2)+w(k)*s*(psi(1,2)*eta(2))*dF;
    C(5,2) = C(5,2)+w(k)*s*(-psi(2,3)*eta(2))*dF;
    C(6,2) = C(6,2)+w(k)*s*(psi(1,3)*eta(2))*dF;
    C(7,2) = C(7,2)+w(k)*(-psi(2,4)*eta(2))*dF;
    C(8,2) = C(8,2)+w(k)*(psi(1,4)*eta(2))*dF;
    %
    C(1,3) = C(1,3)+w(k)*s*(-psi(2,1)*eta(3))*dF;
    C(2,3) = C(2,3)+w(k)*s*(psi(1,1)*eta(3))*dF;
    C(3,3) = C(3,3)+w(k)*s*(-psi(2,2)*eta(3))*dF;
    C(4,3) = C(4,3)+w(k)*s*(psi(1,2)*eta(3))*dF;
    C(5,3) = C(5,3)+w(k)*s*(-psi(2,3)*eta(3))*dF;
    C(6,3) = C(6,3)+w(k)*s*(psi(1,3)*eta(3))*dF;
    C(7,3) = C(7,3)+w(k)*(-psi(2,4)*eta(3))*dF;
    C(8,3) = C(8,3)+w(k)*(psi(1,4)*eta(3))*dF;
end
% -- Matrix u*div(Tau)
B=sparse(8,2);
%
B(1,1) = 2*AT*s/dF;
B(2,2) = 2*AT*s/dF;
B(3,1) = 4/sqrt(2)*AT*s/dF;
B(4,2) = 4/sqrt(2)*AT*s/dF;
B(5,1) = 2*AT*s/dF;
B(6,2) = 2*AT*s/dF;
return