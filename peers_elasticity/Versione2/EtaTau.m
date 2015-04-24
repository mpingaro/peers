function C = EtaTau(DF,DFF,JF,s)

% QUADRATURA DI GAUSS (grado di precisione 5)
a= 9/80;
b=(155+sqrt(15))/2400;
c=(155-sqrt(15))/2400;
w(1)=a; w(2)=b; w(3)=b; w(4)=b;
w(5)=c; w(6)=c; w(7)=c;
q1=1/3; q2=(6+sqrt(15))/21; q3=(9-2*sqrt(15))/21;
q4=(6-sqrt(15))/21; q5=(9+2*sqrt(15))/21;
xnod=[q1 q2 q3 q2 q4 q5 q4]; ynod=[q1 q2 q2 q3 q4 q4 q5];

C = zeros(8,3);
for k=1:7
    x(1) = xnod(k); x(2) = ynod(k);
    % Moltiplicatori recupero simmetria
    eta(1)= 1-x(1)-x(2);
    eta(2)= x(1);
    eta(3)= x(2);
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
    
    % MATRIX C c( eta,as(tau) )
    C(1,1) = C(1,1)+w(k)*s*(-psi(2,1)*eta(1));
    C(2,1) = C(2,1)+w(k)*s*(psi(1,1)*eta(1));
    C(3,1) = C(3,1)+w(k)*s*(-psi(2,2)*eta(1));
    C(4,1) = C(4,1)+w(k)*s*(psi(1,2)*eta(1));
    C(5,1) = C(5,1)+w(k)*s*(-psi(2,3)*eta(1));
    C(6,1) = C(6,1)+w(k)*s*(psi(1,3)*eta(1));
    C(7,1) = C(7,1)+w(k)*(-psi(2,4)*eta(1))*JF;
    C(8,1) = C(8,1)+w(k)*(psi(1,4)*eta(1))*JF;
    %  
    C(1,2) = C(1,2)+w(k)*s*(-psi(2,1)*eta(2));
    C(2,2) = C(2,2)+w(k)*s*(psi(1,1)*eta(2));
    C(3,2) = C(3,2)+w(k)*s*(-psi(2,2)*eta(2));
    C(4,2) = C(4,2)+w(k)*s*(psi(1,2)*eta(2));
    C(5,2) = C(5,2)+w(k)*s*(-psi(2,3)*eta(2));
    C(6,2) = C(6,2)+w(k)*s*(psi(1,3)*eta(2));
    C(7,2) = C(7,2)+w(k)*(-psi(2,4)*eta(2))*JF;
    C(8,2) = C(8,2)+w(k)*(psi(1,4)*eta(2))*JF;
    %
    C(1,3) = C(1,3)+w(k)*s*(-psi(2,1)*eta(3));
    C(2,3) = C(2,3)+w(k)*s*(psi(1,1)*eta(3));
    C(3,3) = C(3,3)+w(k)*s*(-psi(2,2)*eta(3));
    C(4,3) = C(4,3)+w(k)*s*(psi(1,2)*eta(3));
    C(5,3) = C(5,3)+w(k)*s*(-psi(2,3)*eta(3));
    C(6,3) = C(6,3)+w(k)*s*(psi(1,3)*eta(3));
    C(7,3) = C(7,3)+w(k)*(-psi(2,4)*eta(3))*JF;
    C(8,3) = C(8,3)+w(k)*(psi(1,4)*eta(3))*JF;
end

return