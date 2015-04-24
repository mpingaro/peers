% ----------------------------------------------------------------------- %
% --------------- MINI ELEMENT FOR INCOMPRESSIBLE MATERIALS --------------%
% ----------------------------------------------------------------------- %
% ------- By Dott. Ing. Marco Pingaro (P.h.D. Students)           ------- %                                   
% ------- University of Pavia                                     ------- %
% -------                                                         ------- %   
% ----------------------------------------------------------------------- %
%
% ------- Displacemente : P1 + B                                  ------- %
% ------- Pressure      : P1                                      ------- %
% ------------------------------------------------------------------------%
clear all; close all; clc;
% INPUT
length = input('Lunghezza della trave = ');
height = input('Altezza della trave = ');
ndx = input('Numero suddivisioni lungo x = ');
ndy = input('Numero suddivisioni lungo y = ');
f(1,1) = input('Inserire forzante di superficie direzione x = '); 
f(2,1) = input('Inserire forzante di superficie direzione y = ');
young = input('Inserire modulo di Young = ');
poisson = 0.5;      % Impostato 
%-------------------------------------------------------------------------%
dx = length/ndx;
dy = height/ndy;
% Carico geometria dati 
[coordinates,nnod]=CoordinatesType2(ndx,ndy,dx,dy);
[element,nelem]=ElementType2(ndx,ndy);
[mc,ngdlu]=CorrispoMC(element,nelem,nnod);
%
ngdlp = nnod;
ngdlt = ngdlu+ngdlp;
mu2 = young/(1+poisson);
%-------------------------------------------------------------------------%
% GLOBAL STIFF MATRIX
A = sparse(ngdlu,ngdlu);
B = sparse(ngdlu,ngdlp);
b = sparse(ngdlt,1);
% ASSEMBLY GLOBAL MATRIX A, B and b (Global Load Vector)  
for k = 1:nelem
    P([1 2])=coordinates(element(k,1),[1 2]);
    P([3 4])=coordinates(element(k,2),[1 2]);
    P([5 6])=coordinates(element(k,3),[1 2]);
    [AELEM,BELEM,JF] = GradGrad(P);
    [bf] = BodyLoad(JF,f);
    for i=1:8
        for j=1:8
           A(mc(k,i),mc(k,j)) = A(mc(k,i),mc(k,j))+mu2*AELEM(i,j);
        end
        b(mc(k,i),1) = b(mc(k,i),1) + bf(i,1);
        for ii = 1:3
            B(mc(k,i),element(k,ii)) = B(mc(k,i),element(k,ii))...
                +BELEM(i,ii);
        end
    end
end
K = [A, -B; B', sparse(nnod,nnod)];
clear A B bf AELEM BELEM JF P
clear young poisson f
%-------------------------------------------------------------------------%
% Vincoli EX : 
%      Trave incernierata e incarrellata ai due estremi bassi
%-------------------------------------------------------------------------%
I = [1, 2, 2*(ndx+1)];
for jj = 1:size(I,2)
     K(I(jj),:) = 0;
     K(I(jj),I(jj)) = 1;
     b(I(jj),1) = 0;
end
% Solution
soluz = K\b;
Ux = soluz(1:2:2*nnod); 
Uy = soluz(2:2:2*nnod);
pressure = soluz(ngdlu+1:ngdlt);
defo =coordinates + [Ux, Uy];
% PLOT SOLUTION
%
figure
plotmesh(element,coordinates)
title('Undeformed Mesh','fontsize',14);
figure
plotmesh(element,defo)
title('Deformed Mesh','fontsize',14);
%
figure
ShowContDisp(element,coordinates,Ux)
title('Displacement Ux','fontsize',14);
%
figure
ShowContDisp(element,coordinates,Uy)
title('Displacement Uy','fontsize',14);
%
figure
ShowContDisp(element,coordinates,pressure)
title('Pressure','fontsize',14);