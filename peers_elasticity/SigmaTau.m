function A = SigmaTau(DF,DFF,JF,mu,lambda,s)
%   mu, lambda = costanti di Lamè
%   s = segno delle funzioni di forma (inversione flussi)

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
for k=1:nquad
    x(1) = xnod(k); x(2) = ynod(k);
    % RT0
    psi(1,1) = DF(1,1)*x(1)+DF(2,1)*(-1+x(2));
    psi(2,1) = DF(3,1)*x(1)+DF(4,1)*(-1+x(2));
    psi(1,2) = DF(1,1)*2/sqrt(2)*x(1)+DF(2,1)*2/sqrt(2)*x(2);
    psi(2,2) = DF(3,1)*2/sqrt(2)*x(1)+DF(4,1)*2/sqrt(2)*x(2);
    psi(1,3) = DF(1,1)*(-1+x(1))+DF(2,1)*x(2);
    psi(2,3) = DF(3,1)*(-1+x(1))+DF(4,1)*x(2);
    % BOLLA
    dBx = (x(2)-x(2)*x(2)-2*x(1)*x(2))*120; 
    dBy = (x(1)-x(1)*x(1)-2*x(1)*x(2))*120;
    dBT(1,1) = DFF(1,1)*dBx+DFF(3,1)*dBy;
    dBT(2,1) = DFF(2,1)*dBx+DFF(4,1)*dBy;
    psi(1,4) =  dBT(2,1)*JF; 
    psi(2,4) = -dBT(1,1)*JF;
    
    % MATRIX A a( C*sigma,tau )
    A(1,1) = A(1,1)+w(k)*(t1*(psi(1,1)*psi(1,1)...
        +psi(2,1)*psi(2,1))-t2*psi(1,1)*psi(1,1))/JF;        
    A(1,2) = A(1,2)+w(k)*(-t2*psi(1,1)*psi(2,1))/JF;
    A(1,3) = A(1,3)+w(k)*(t1*(psi(1,1)*psi(1,2)...
        +psi(2,1)*psi(2,2))-t2*psi(1,1)*psi(1,2))/JF;
    A(1,4) = A(1,4)+w(k)*(-t2*psi(1,1)*psi(2,2))/JF;
    A(1,5) = A(1,5)+w(k)*(t1*(psi(1,1)*psi(1,3)...
        +psi(2,1)*psi(2,3))-t2*psi(1,1)*psi(1,3))/JF;
    A(1,6) = A(1,6)+w(k)*(-t2*psi(1,1)*psi(2,3))/JF;
    %
    A(2,2) = A(2,2)+w(k)*(t1*(psi(1,1)*psi(1,1)...
        +psi(2,1)*psi(2,1))-t2*psi(2,1)*psi(2,1))/JF;
    A(2,3) = A(2,3)+w(k)*(-t2*psi(2,1)*psi(1,2))/JF;
    A(2,4) = A(2,4)+w(k)*(t1*(psi(1,1)*psi(1,2)...
        +psi(2,1)*psi(2,2))-t2*psi(2,1)*psi(2,2))/JF;
    A(2,5) = A(2,5)+w(k)*(-t2*psi(2,1)*psi(1,3))/JF;
    A(2,6) = A(2,6)+w(k)*(t1*(psi(1,1)*psi(1,3)...
        +psi(2,1)*psi(2,3))-t2*psi(2,1)*psi(2,3))/JF;
    %
    A(3,3) = A(3,3)+w(k)*(t1*(psi(1,2)*psi(1,2)...
        +psi(2,2)*psi(2,2))-t2*psi(1,2)*psi(1,2))/JF;
    A(3,4) = A(3,4)+w(k)*(-t2*psi(1,2)*psi(2,2))/JF;
    A(3,5) = A(3,5)+w(k)*(t1*(psi(1,2)*psi(1,3)...
        +psi(2,2)*psi(2,3))-t2*psi(1,2)*psi(2,3))/JF;
    A(3,6) = A(3,6)+w(k)*(-t2*psi(1,2)*psi(2,3))/JF;
    %
    A(4,4) = A(4,4)+w(k)*(t1*(psi(1,2)*psi(1,2)...
        +psi(2,2)*psi(2,2))-t2*psi(2,2)*psi(2,2))/JF;
    A(4,5) = A(4,5)+w(k)*(-t2*psi(2,2)*psi(1,3))/JF;
    A(4,6) = A(4,6)+w(k)*(t1*(psi(1,2)*psi(1,3)...
        +psi(2,2)*psi(2,3))-t2*psi(2,2)*psi(2,3))/JF;
    %
    A(5,5) = A(5,5)+w(k)*(t1*(psi(1,3)*psi(1,3)...
        +psi(2,3)*psi(2,3))-t2*psi(1,3)*psi(1,3))/JF;
    A(5,6) = A(5,6)+w(k)*(-t2*psi(1,3)*psi(2,3))/JF;
    
    A(6,6) = A(6,6)+w(k)*(t1*(psi(1,3)*psi(1,3)...
        +psi(2,3)*psi(2,3))-t2*psi(2,3)*psi(2,3))/JF;
    %
    A(7,7) = A(7,7)+w(k)*(t1*(psi(1,4)*psi(1,4)...
        +psi(2,4)*psi(2,4))-t2*psi(1,4)*psi(1,4))*JF;
    A(7,8) = A(7,8)+w(k)*(-t2*psi(1,4)*psi(2,4))*JF;
    A(8,8) = A(8,8)+w(k)*(t1*(psi(1,4)*psi(1,4)...
        +psi(2,4)*psi(2,4))-t2*psi(2,4)*psi(2,4))*JF;
    %
    A(1,7) = A(1,7)+w(k)*s*(t1*(psi(1,1)*psi(1,4)... 
        +psi(2,1)*psi(2,4))-t2*psi(1,1)*psi(1,4)); 
    A(1,8) = A(1,8)+w(k)*s*(-t2*psi(1,1)*psi(2,4));
    A(2,7) = A(2,7)+w(k)*s*(-t2*psi(2,1)*psi(1,4));
    A(2,8) = A(2,8)+w(k)*s*(t1*(psi(1,1)*psi(1,4)...
        +psi(2,1)*psi(2,4))-t2*psi(2,1)*psi(2,4));
    A(3,7) = A(3,7)+w(k)*s*(t1*(psi(1,2)*psi(1,4)...
        +psi(2,2)*psi(2,4))-t2*psi(1,2)*psi(1,4));
    A(3,8) = A(3,8)+w(k)*s*(-t2*psi(1,2)*psi(2,4)); 
    A(4,7) = A(4,7)+w(k)*s*(-t2*psi(2,2)*psi(1,4));
    A(4,8) = A(4,8)+w(k)*s*(t1*(psi(1,2)*psi(1,4)...
        +psi(2,2)*psi(2,4))-t2*psi(2,2)*psi(2,4));
    A(5,7) = A(5,7)+w(k)*s*(t1*(psi(1,3)*psi(1,4)...
        +psi(2,3)*psi(2,4))-t2*psi(1,3)*psi(1,4));
    A(5,8) = A(5,8)+w(k)*s*(-t2*psi(1,3)*psi(2,4));
    A(6,7) = A(6,7)+w(k)*s*(-t2*psi(2,3)*psi(1,4));
    A(6,8) = A(6,8)+w(k)*s*(t1*(psi(1,3)*psi(1,4)...
        +psi(2,3)*psi(2,4))-t2*psi(2,3)*psi(2,4));
end
return
