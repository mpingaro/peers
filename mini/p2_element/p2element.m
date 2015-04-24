function A= p2element(P,nu,E)

%% QUADRATURA DI GAUSS (grado di precisione 5)
% Weight of quadrature
w(1,1)=0.112500000000000; w(1,2)=0.066197076394253; w(1,3)=0.066197076394253; 
w(1,4)=0.066197076394253; w(1,5)=0.062969590272414; w(1,6)=0.062969590272414; 
w(1,7)=0.062969590272414;
% Point value 
gauss_x(1,1) = 0.333333333333333; gauss_x(1,2) = 0.470142064105115; 
gauss_x(1,3) = 0.059715871789770; gauss_x(1,4) = 0.470142064105115; 
gauss_x(1,5) = 0.101286507323456; gauss_x(1,6) = 0.797426985353087;
gauss_x(1,7) = 0.101286507323456; 
gauss_y(1,1) = 0.333333333333333; gauss_y(1,2) = 0.470142064105115; 
gauss_y(1,3) = 0.470142064105115; gauss_y(1,4) = 0.059715871789770;
gauss_y(1,5) = 0.101286507323456; gauss_y(1,6) = 0.101286507323456;
gauss_y(1,7) = 0.797426985353087;
nqd = 7;
%% - Legame    -----------------------------------------------------------%
lambda = E*nu/( (1+nu)*(1-2*nu) ); 
G = E/(2*(1+nu));
C = [lambda+2*G, lambda, 0;lambda, lambda+2*G, 0; 0, 0, G];

%% STIFFNESS MATRIX
A = sparse(12,12);
for k =1:nqd
    x = gauss_x(1,k); y = gauss_y(1,k);
     
    DF(1,1) = (4*x+4*y-3)*P(1,1) + (4*x-1)*P(1,3) + (4-8*x-4*y)*P(1,7)...
        +4*y*P(1,9) - 4*y*P(1,11);
    DF(1,2) = (4*x+4*y-3)*P(1,1) + (4*y-1)*P(1,5) -4*x*P(1,7)...
        +4*x*P(1,9) + (4-4*x-8*y)*P(1,11);
    DF(2,1) = (4*x+4*y-3)*P(1,2) + (4*x-1)*P(1,4) + (4-8*x-4*y)*P(1,8)...
        +4*y*P(1,10) - 4*y*P(1,12);
    DF(2,2) = (4*x+4*y-3)*P(1,2) + (4*y-1)*P(1,6) -4*x*P(1,8)...
        +4*x*P(1,10) + (4-4*x-8*y)*P(1,12);
    
    JF = DF(1,1)*DF(2,2)-DF(1,2)*DF(2,1);
    %% Inverse of Jacobian Transformation
    DFF(1,1) = DF(2,2)/JF;
    DFF(1,2) = -DF(2,1)/JF;
    DFF(2,1) = -DF(1,1)/JF;
    DFF(2,2) = DF(1,1)/JF;
    %% Grad Shape Function in reference domain
    grad(1,1) = 4*x+4*y-3; grad(2,1) = 4*x+4*y-3;
    grad(1,2) = 4*x-1;     grad(2,2) = 0;
    grad(1,3) = 0;         grad(2,3) = 4*y-1;
    grad(1,4) = 4-8*x-4*y; grad(2,4) = -4*x;
    grad(1,5) = 4*y;       grad(2,5) = 4*x;
    grad(1,6) = -4*y;      grad(2,6) = 4-4*x-8*y;
    %% Grad Shape Function in Physical Domain
    gradT(:,1) = DFF*grad(:,1); 
    gradT(:,2) = DFF*grad(:,2);
    gradT(:,3) = DFF*grad(:,3); 
    gradT(:,4) = DFF*grad(:,4);
    gradT(:,5) = DFF*grad(:,5); 
    gradT(:,6) = DFF*grad(:,6);
    %% Shepe Function (epsiXX,epsiYY,gammaXY)
    epsi(:,1)  = [gradT(1,1); 0; gradT(2,1)]; 
    epsi(:,2)  = [0; gradT(2,1); gradT(1,1)];
    epsi(:,3)  = [gradT(1,2); 0; gradT(2,2)];
    epsi(:,4)  = [0; gradT(2,2); gradT(1,2)];
    epsi(:,5)  = [gradT(1,3); 0; gradT(2,3)];
    epsi(:,6)  = [0; gradT(2,3); gradT(1,3)];
    epsi(:,7)  = [gradT(1,4); 0; gradT(2,4)]; 
    epsi(:,8)  = [0; gradT(2,4); gradT(1,4)];
    epsi(:,9)  = [gradT(1,5); 0; gradT(2,5)];
    epsi(:,10) = [0; gradT(2,5); gradT(1,5)];
    epsi(:,11) = [gradT(1,6); 0; gradT(2,6)];
    epsi(:,12) = [0; gradT(2,6); gradT(1,6)];
    
    for i=1:12
        for j=1:12
            sigma = C*epsi(:,i);
            A(i,j) = A(i,j) + w(1,k)*(sigma(1,1)*epsi(1,j) +...
                sigma(2,1)*epsi(2,j) + sigma(3,1)*epsi(3,j))*JF;
        end
    end
      
end

return