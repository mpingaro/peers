% ----------------------------------------------------------------------- %
%--  Matrix of point value (parent triangle)                            --%
%--                                                                     --%
%--                                                                     --%
%-- a1 + a2*x + a3*y + a4*x^2 + a5*x*y + a6*y^2 --> Shape Function      --%
%--                                                                     --%
%-- coeff(6,6) matrix a(:,i)                                            --%
% ----------------------------------------------------------------------- %
clear all;close all; clc;
syms x y


P = [0, 0; 1, 0; 0, 1; 0.5, 0; 0.5, 0.5; 0, 0.5];

coeff = zeros(6,6);
for i = 1:6
    K = zeros(6,6);
    for j = 1:6
        K(j,:) = [1, P(j,1), P(j,2), P(j,1)*P(j,1),...
            P(j,1)*P(j,2), P(j,2)*P(j,2)];
    end
    b(:,1) = zeros(6,1);
    b(i,1) = 1;
    coeff(:,i) = K\b;
end

psi(1,1) = coeff(:,1)'*[1; x; y; x.^2; x.*y; y.^2];
psi(2,1) = coeff(:,2)'*[1; x; y; x.^2; x.*y; y.^2];
psi(3,1) = coeff(:,3)'*[1; x; y; x.^2; x.*y; y.^2];
psi(4,1) = coeff(:,4)'*[1; x; y; x.^2; x.*y; y.^2];
psi(5,1) = coeff(:,5)'*[1; x; y; x.^2; x.*y; y.^2];
psi(6,1) = coeff(:,6)'*[1; x; y; x.^2; x.*y; y.^2];

disp(psi);