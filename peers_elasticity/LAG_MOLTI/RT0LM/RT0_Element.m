function [BELEM,CELEM] = RT0_Element(msh,point)

%% MATRIX C
CELEM = zeros(3,1);
CELEM(3,1) = 2*msh.area;

%% MATRIX B
BELEM = zeros(3,3);

BELEM(1,1) = msh.area; 
BELEM(2,2) = msh.area;

%% QUADRATURA DI GAUSS (grado di precisione 5)
% Weight of quadrature
weight(1,1) = 0.112500000000000; weight(1,2) = 0.066197076394253; 
weight(1,3) = 0.066197076394253; weight(1,4) = 0.066197076394253; 
weight(1,5) = 0.062969590272414; weight(1,6) = 0.062969590272414; 
weight(1,7) = 0.062969590272414;
% Point value 
gauss(1,1) = 0.333333333333333; gauss(1,2) = 0.470142064105115; 
gauss(1,3) = 0.059715871789770; gauss(1,4) = 0.470142064105115; 
gauss(1,5) = 0.101286507323456; gauss(1,6) = 0.797426985353087;
gauss(1,7) = 0.101286507323456; 

gauss(2,1) = 0.333333333333333; gauss(2,2) = 0.470142064105115; 
gauss(2,3) = 0.470142064105115; gauss(2,4) = 0.059715871789770;
gauss(2,5) = 0.101286507323456; gauss(2,6) = 0.101286507323456;
gauss(2,7) = 0.797426985353087;
nqd = 7;

gauss = msh.jac*gauss + [point(1,1)*ones(1,7);point(1,2)*ones(1,7)];

for k=1:nqd
    x(1) = gauss(1,k); x(2) = gauss(2,k);
    
    % Shape
    eta(1) = x(1)-msh.bari(1);
    eta(2) = x(2)-msh.bari(2);

    BELEM(3,3) = BELEM(3,3) + weight(k)*( eta(1)^2 + eta(2)^2 )*msh.djac;

end
% Versione Carstersen
% s = norm(point(2,:)-point(1,:))^2 + norm(point(3,:)-point(2,:))^2 + norm(point(3,:)-point(1,:))^2;
% BB = 1/36*msh.area*s;

end