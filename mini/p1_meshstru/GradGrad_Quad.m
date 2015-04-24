function [AELEM,bf] = GradGrad_Quad(point,E,nu,f)  

%% Gauss Quadrature
%  4 Point Quadrature
%  Weight of Quadrature is uqual to 1

% Direction x
gauss_x(1,1) = -0.577350269189626; gauss_x(1,2) =  0.577350269189626;
gauss_x(1,3) =  0.577350269189626; gauss_x(1,4) = -0.577350269189626;
% Direction y
gauss_y(1,1) = -0.577350269189626; gauss_y(1,2) = -0.577350269189626;
gauss_y(1,3) =  0.577350269189626; gauss_y(1,4) =  0.577350269189626;

npg = 4;

%% ELEMENTARY MATRIX A & B
AELEM = zeros(8,8); bf = zeros(8,1);
for k = 1:npg
    x = gauss_x(1,k); y = gauss_y(1,k);
    w = 1;
    
    J(1,1) = 0.25*(1-y)*(point(2,1)-point(1,1))...
        + 0.25*(1+y)*(point(3,1)-point(4,1));       % x_u
    J(1,2) = 0.25*(1-x)*(point(4,1)-point(1,1))...
        + 0.25*(1+x)*(point(3,1)-point(2,1));       % x_v
    J(2,1) = 0.25*(1-y)*(point(2,2)-point(1,2))...
        + 0.25*(1+y)*(point(3,2)-point(4,2));       % y_u
    J(2,2) = 0.25*(1-x)*(point(4,2)-point(1,2))...
        + 0.25*(1+x)*(point(3,2)-point(2,2));       % y_v
    
    DJ = J(1,1)*J(2,2)-J(1,2)*J(2,1);               % Det (J)
    
    % Inverse and transpose of Jacobian Matrix
    JJ(1,1) = J(2,2)/DJ;  
    JJ(1,2) = -J(2,1)/DJ;
    JJ(2,1) = -J(1,2)/DJ; 
    JJ(2,2) = J(1,1)/DJ;
    % --------------------------------------------------------------------%
    % epsiXX, epsiYY, gammaXY    
    grad(:,1) = JJ*[-0.25*(1-y); -0.25*(1-x)];
    grad(:,2) = JJ*[0.25*(1-y); -0.25*(1+x)];
    grad(:,3) = JJ*[0.25*(1+y); 0.25*(1+x)];
    grad(:,4) = JJ*[-0.25*(1+y); 0.25*(1-x)];
    %
    grad = JJ*grad;
    
    % - Legame    --------------------------------------------------------%
    lambda = E*nu/( (1+nu)*(1-2*nu) ); 
    G = E/(2*(1+nu));
    C = [lambda+2*G, lambda, 0;lambda, lambda+2*G, 0; 0, 0, G];
    
    % - Shape functions 1-8    -------------------------------------------%
    N =[psi(1) 0; 
        0 psi(1); 
        psi(2) 0;
        0 psi(2); 
        psi(3) 0; 
        0 psi(3); 
        psi(4) 0; 
        0 psi(4)]; 
    
    B = [grad(1,1) 0 grad(1,2) 0 grad(1,3) 0 grad(1,4) 0;
        0 grad(2,1) 0 grad(2,2) 0 grad(2,3) 0 grad(2,4);
        grad(2,1) grad(1,1) grad(2,2) grad(1,2) grad(2,3) grad(1,3) grad(2,4) grad(1,4)];
    
    AELEM = AELEM + w .* (B' * C * B) .* DJ;     % Stiffness matrix
    bf(:,1) = bf(:,1) + w.* (N*f) .*DJ;          % Bulk load
    
end

end