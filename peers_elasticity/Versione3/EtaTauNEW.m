function C = EtaTauNEW(psiF,etaF,JF,s)

% QUADRATURA DI GAUSS (grado di precisione 5)
a= 9/80;
b=(155+sqrt(15))/2400;
c=(155-sqrt(15))/2400;
w(1,1)=a; w(1,2)=b; w(1,3)=b; w(1,4)=b;
w(1,5)=c; w(1,6)=c; w(1,7)=c;

C = sparse(8,3);
% MATRIX C c( eta,as(tau) )
C(1,1) = s*(-w(1,:).*psiF(2,:)*etaF(1,:)');
C(2,1) = s*( w(1,:).*psiF(1,:)*etaF(1,:)');
C(3,1) = s*(-w(1,:).*psiF(4,:)*etaF(1,:)');
C(4,1) = s*( w(1,:).*psiF(3,:)*etaF(1,:)');
C(5,1) = s*(-w(1,:).*psiF(6,:)*etaF(1,:)');
C(6,1) = s*( w(1,:).*psiF(5,:)*etaF(1,:)');
C(7,1) = s*(-w(1,:).*psiF(8,:)*etaF(1,:)')*JF;
C(8,1) = s*( w(1,:).*psiF(7,:)*etaF(1,:)')*JF;
%  
C(1,2) = s*(-w(1,:).*psiF(2,:)*etaF(2,:)');
C(2,2) = s*( w(1,:).*psiF(1,:)*etaF(2,:)');
C(3,2) = s*(-w(1,:).*psiF(4,:)*etaF(2,:)');
C(4,2) = s*( w(1,:).*psiF(3,:)*etaF(2,:)');
C(5,2) = s*(-w(1,:).*psiF(6,:)*etaF(2,:)');
C(6,2) = s*( w(1,:).*psiF(5,:)*etaF(2,:)');
C(7,2) = s*(-w(1,:).*psiF(8,:)*etaF(2,:)')*JF;
C(8,2) = s*( w(1,:).*psiF(7,:)*etaF(2,:)')*JF;
%
C(1,3) = s*(-w(1,:).*psiF(2,:)*etaF(3,:)');
C(2,3) = s*( w(1,:).*psiF(1,:)*etaF(3,:)');
C(3,3) = s*(-w(1,:).*psiF(4,:)*etaF(3,:)');
C(4,3) = s*( w(1,:).*psiF(3,:)*etaF(3,:)');
C(5,3) = s*(-w(1,:).*psiF(6,:)*etaF(3,:)');
C(6,3) = s*( w(1,:).*psiF(5,:)*etaF(3,:)');
C(7,3) = s*(-w(1,:).*psiF(8,:)*etaF(3,:)')*JF;
C(8,3) = s*( w(1,:).*psiF(7,:)*etaF(3,:)')*JF;

return