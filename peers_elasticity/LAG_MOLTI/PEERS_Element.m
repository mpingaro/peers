    function [BELEM,CELEM,DELEM] = PEERS_Element(msh,point,alpha,beta)

%% MATRIX C --------------------------------------------------------------%
CELEM = zeros(8,2);

CELEM(5,1) = 2*msh.area;
CELEM(5,2) = CELEM(5,1);
CELEM(6,1) = CELEM(5,1);
CELEM(6,2) = -CELEM(5,1);

%% MATRIX D --------------------------------------------------------------%
DELEM = zeros(8,3);

DELEM(4,1) = -2*msh.area/3; 
DELEM(4,2) = DELEM(4,1);
DELEM(4,3) = DELEM(4,1);
%
for i=1:3
    DELEM(5,i) = msh.area*( point(i,2)-msh.bari(2)-point(i,1)+msh.bari(1) )/12;
    DELEM(6,i) = msh.area*( point(i,2)-msh.bari(2)+point(i,1)-msh.bari(1) )/12;
end
%
DELEM(7,1) = ( point(2,2)-point(3,2)+point(3,1)-point(2,1) )/120;
DELEM(7,2) = ( point(3,2)-point(1,2)+point(1,1)-point(3,1) )/120;
DELEM(7,3) = ( point(1,2)-point(2,2)+point(2,1)-point(1,1) )/120;
%
DELEM(8,1) = ( point(2,2)-point(3,2)-point(3,1)+point(2,1) )/120;
DELEM(8,2) = ( point(3,2)-point(1,2)-point(1,1)+point(3,1) )/120;
DELEM(8,3) = ( point(1,2)-point(2,2)-point(2,1)+point(1,1) )/120;

%% MATRIX B --------------------------------------------------------------%
BELEM = zeros(8,8);

BELEM(1,1) = 2*msh.area*(alpha + 2*beta); 
BELEM(2,2) = 2*msh.area*alpha;
BELEM(3,3) = BELEM(2,2);
BELEM(4,4) = BELEM(2,2);

BELEM(5,8) = -msh.area*beta/30; 
BELEM(8,5) = BELEM(5,8);
BELEM(6,7) = -BELEM(5,8);
BELEM(7,6) = BELEM(6,7);

[gauss, weight] = Gauss_Quadrature(2);
gauss = msh.jac*gauss + [point(1,1)*ones(1,3);point(1,2)*ones(1,3)];
for k=1:3
    x(1) = gauss(1,k); x(2) = gauss(2,k);
    % Shape
    eta(1) = x(1) - msh.bari(1);
    eta(2) = x(2) - msh.bari(2);
    
    BELEM(5,5) = BELEM(5,5) + weight(k)*( 2*alpha*( eta(1)^2 + eta(2)^2 )...
       - beta*( eta(1)+eta(2) )^2 );%*msh.djac;
    BELEM(5,6) = BELEM(5,6) + weight(k)*( -beta*( eta(1)+eta(2) )...
       *( eta(1)-eta(2) ) );%*msh.djac;
    BELEM(6,6) = BELEM(6,6) + weight(k)*( 2*alpha*( eta(1)^2 + eta(2)^2 )...
       - beta*( eta(1)-eta(2) )^2 );%*msh.djac;
end
    
[gauss, weight] = Gauss_Quadrature(5);
%gauss = msh.jac*gauss + [point(1,1)*ones(1,7);point(1,2)*ones(1,7)];
for k=1:7
    x(1) = gauss(1,k); x(2) = gauss(2,k); 
    % BOLLA
    dBx = (x(2)-x(2)*x(2)-2*x(1)*x(2))/120; 
    dBy = (x(1)-x(1)*x(1)-2*x(1)*x(2))/120;
    % Trasformo il gradiente dal parametrico al fisico.
    grad_B = inv(msh.jac)'*[dBx;dBy];
    curl_B = [grad_B(2); -grad_B(1)];
    
    BELEM(7,7) = BELEM(7,7) + weight(k)*( 2*alpha*( curl_B(1)^2 + curl_B(2)^2 )...
        - beta*( curl_B(1)+curl_B(2) )^2 )*msh.djac;
    BELEM(7,8) = BELEM(7,8) + weight(k)*( - beta*( curl_B(1)+curl_B(2) )...
        *( curl_B(1)-curl_B(2) ) )*msh.djac;
    BELEM(8,8) = BELEM(8,8) + weight(k)*( 2*alpha*( curl_B(1)^2 + curl_B(2)^2 )...
        - beta*( curl_B(1)-curl_B(2) )^2 )*msh.djac;

end
BELEM(6,5) = BELEM(5,6);
BELEM(8,7) = BELEM(7,8);

% CALCOLO CARSTENSEN
% BELEM66 = ( (2*alpha+beta)/12*( norm(point(1,:)-msh.bari)^2 +...
%     norm(point(2,:)-msh.bari)^2 + norm(point(3,:)-msh.bari)^2 )...
%     -beta/6*( (point(1,1)-msh.bari(1) )*( point(1,2)-msh.bari(2) ) +...
%     (point(2,1)-msh.bari(1) )*( point(2,2)-msh.bari(2) )+...
%     (point(3,1)-msh.bari(1) )*( point(3,2)-msh.bari(2) ) ) )*msh.area;
%     

% BELEM56 = (beta/12)*( (point(1,1)-msh.bari(1) )^2 - ( point(1,2)-msh.bari(2) )^2+...
%     ( point(2,1)-msh.bari(1) )^2 - ( point(2,2)-msh.bari(2) )^2+...
%     ( point(3,1)-msh.bari(1) )^2 - ( point(3,2)-msh.bari(2) )^2 )*msh.area;

% disp(BELEM66);
% disp(BELEM(6,6));
% disp(BELEM56);
% disp( BELEM(5,6) );

end