function [psiF,etaF] = EvaluateShapePhysics(DF,DFF,JF,psi,eta)

etaF = eta;

for i=1:7
   psiF(1,i) =  DF(1,1)*psi(1,i)+DF(2,1)*psi(2,i);
   psiF(2,i) =  DF(3,1)*psi(1,i)+DF(4,1)*psi(2,i);
   psiF(3,i) =  DF(1,1)*psi(3,i)+DF(2,1)*psi(4,i);
   psiF(4,i) =  DF(3,1)*psi(3,i)+DF(4,1)*psi(4,i);
   psiF(5,i) =  DF(1,1)*psi(5,i)+DF(2,1)*psi(6,i);
   psiF(6,i) =  DF(3,1)*psi(5,i)+DF(4,1)*psi(6,i);
   psiF(7,i) =  DFF(2,1)*psi(7,i)*JF+DFF(4,1)*psi(8,i)*JF;
   psiF(8,i) = -DFF(1,1)*psi(7,i)*JF-DFF(3,1)*psi(8,7)*JF; 
end
return