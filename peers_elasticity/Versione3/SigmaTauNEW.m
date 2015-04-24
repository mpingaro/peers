function A = SigmaTauNEW(psiF,JF,mu,lambda,s)

% QUADRATURA DI GAUSS (grado di precisione 5)
a= 9/80;
b=(155+sqrt(15))/2400;
c=(155-sqrt(15))/2400;
w(1,1)=a; w(1,2)=b; w(1,3)=b; w(1,4)=b;
w(1,5)=c; w(1,6)=c; w(1,7)=c;
%
% Elastic constant
t1 = 1/(2*mu); t2 = lambda/(4*mu*(mu+lambda));
A = sparse(8,8);                         % Storo solo triangolare superiore
% MATRIX A a( C*sigma,tau )
A(1,1) = t1*(w(1,:).*psiF(1,:)*psiF(1,:)'+w(1,:).*psiF(2,:)*psiF(2,:)')/JF...
    -t2*(w(1,:).*psiF(1,:)*psiF(1,:)')/JF;
A(1,2) = -t2*w(1,:).*psiF(1,:)*psiF(2,:)'/JF; 
A(1,3) = t1*(w(1,:).*psiF(1,:)*psiF(3,:)'+w(1,:).*psiF(2,:)*psiF(4,:)')/JF...
    -t2*(w(1,:).*psiF(1,:)*psiF(3,:)')/JF;
A(1,4) = -t2*(w(1,:).*psiF(1,:)*psiF(4,:)')/JF;
A(1,5) = t1*(w(1,:).*psiF(1,:)*psiF(5,:)'+w(1,:).*psiF(2,:)*psiF(6,:)')/JF...
    -t2*(w(1,:).*psiF(1,:)*psiF(5,:)')/JF;
A(1,6) = -t2*(w(1,:).*psiF(1,:)*psiF(6,:)')/JF;
A(1,7) = t1*s*(w(1,:).*psiF(1,:)*psiF(7,:)'+w(1,:).*psiF(2,:)*psiF(8,:)')...
    -t2*s*(w(1,:).*psiF(1,:)*psiF(7,:)');
A(1,8) = -t2*s*(w(1,:).*psiF(1,:)*psiF(8,:)');
%
A(2,2) = t1*(w(1,:).*psiF(1,:)*psiF(1,:)'+w(1,:).*psiF(2,:)*psiF(2,:)')/JF...
    -t2*(w(1,:).*psiF(2,:)*psiF(2,:)')/JF;
A(2,3) = -t2*(w(1,:).*psiF(2,:)*psiF(3,:)')/JF;
A(2,4) = t1*(w(1,:).*psiF(1,:)*psiF(3,:)'+w(1,:).*psiF(2,:)*psiF(4,:)')/JF...
    -t2*(w(1,:).*psiF(2,:)*psiF(4,:)')/JF;
A(2,5) = -t2*(w(1,:).*psiF(2,:)*psiF(5,:)')/JF;
A(2,6) = t1*(w(1,:).*psiF(1,:)*psiF(5,:)'+w(1,:).*psiF(2,:)*psiF(6,:)')/JF...
    -t2*(w(1,:).*psiF(2,:)*psiF(6,:)')/JF;    
A(2,7) = -t2*s*(w(1,:).*psiF(2,:)*psiF(7,:)');
A(2,8) = t1*s*(w(1,:).*psiF(1,:)*psiF(7,:)'+w(1,:).*psiF(2,:)*psiF(8,:)')...
    -t2*s*(w(1,:).*psiF(2,:)*psiF(8,:)');
%
A(3,3) = t1*(w(1,:).*psiF(3,:)*psiF(3,:)'+w(1,:).*psiF(4,:)*psiF(4,:)')/JF...
    -t2*(w(1,:).*psiF(3,:)*psiF(3,:)')/JF;
A(3,4) = -t2*(w(1,:).*psiF(3,:)*psiF(4,:)')/JF;
A(3,5) = t1*(w(1,:).*psiF(3,:)*psiF(5,:)'+w(1,:).*psiF(4,:)*psiF(6,:)')/JF...
    -t2*(w(1,:).*psiF(3,:)*psiF(5,:)')/JF; % DA CONTROLLARE psiF(6,:) psi(2,3) prima!!!! forse è psi(1,3)
A(3,6) = -t2*(w(1,:).*psiF(3,:)*psiF(6,:)')/JF;
A(3,7) = t1*s*(w(1,:).*psiF(3,:)*psiF(7,:)'+w(1,:).*psiF(4,:)*psiF(8,:)')...
    -t2*s*(w(1,:).*psiF(3,:)*psiF(7,:)');
A(3,8) = -t2*s*(w(1,:).*psiF(3,:)*psiF(8,:)');
%
A(4,4) = t1*(w(1,:).*psiF(3,:)*psiF(3,:)'+w(1,:).*psiF(4,:)*psiF(4,:)')/JF...
    -t2*(w(1,:).*psiF(4,:)*psiF(4,:)')/JF;
A(4,5) = -t2*(w(1,:).*psiF(4,:)*psiF(5,:)')/JF;
A(4,6) = t1*(w(1,:).*psiF(3,:)*psiF(5,:)'+w(1,:).*psiF(4,:)*psiF(6,:)')/JF...
    -t2*(w(1,:).*psiF(4,:)*psiF(6,:)')/JF;
A(4,7) = -t2*s*(w(1,:).*psiF(4,:)*psiF(7,:)');
A(4,8) = t1*s*(w(1,:).*psiF(3,:)*psiF(7,:)'+w(1,:).*psiF(4,:)*psiF(8,:)')...
    -t2*s*(w(1,:).*psiF(4,:)*psiF(8,:)');
%
A(5,5) = t1*(w(1,:).*psiF(5,:)*psiF(5,:)'+w(1,:).*psiF(6,:)*psiF(6,:)')/JF...
    -t2*(w(1,:).*psiF(5,:)*psiF(5,:)')/JF;
A(5,6) = -t2*(w(1,:).*psiF(5,:)*psiF(6,:)')/JF;
A(5,7) = t1*s*(w(1,:).*psiF(5,:)*psiF(7,:)'+w(1,:).*psiF(6,:)*psiF(8,:)')...
    -t2*s*(w(1,:).*psiF(5,:)*psiF(7,:)');
A(5,8) = -t2*s*(w(1,:).*psiF(5,:)*psiF(8,:)');
%
A(6,6) = t1*(w(1,:).*psiF(5,:)*psiF(5,:)'+w(1,:).*psiF(6,:)*psiF(6,:)')/JF...
    -t2*(w(1,:).*psiF(6,:)*psiF(6,:)')/JF;
A(6,7) = -t2*s*(w(1,:).*psiF(6,:)*psiF(7,:)');
A(6,8) = t1*s*(w(1,:).*psiF(5,:)*psiF(7,:)'+w(1,:).*psiF(6,:)*psiF(8,:)')...
    -t2*s*(w(1,:).*psiF(6,:)*psiF(8,:)');
%
A(7,7) = t1*(w(1,:).*psiF(7,:)*psiF(7,:)'+w(1,:).*psiF(8,:)*psiF(8,:)')*JF...
    -t2*(w(1,:).*psiF(7,:)*psiF(7,:)')*JF;    
A(7,8) = -t2*(w(1,:).*psiF(7,:)*psiF(8,:)')*JF;
%
A(8,8) = t1*(w(1,:).*psiF(7,:)*psiF(7,:)'+w(1,:).*psiF(8,:)*psiF(8,:)')*JF...
    -t2*(w(1,:).*psiF(8,:)*psiF(8,:)')*JF;
return