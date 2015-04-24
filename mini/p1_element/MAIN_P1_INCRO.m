% ----------------------------------------------------------------------- %
% ---------------- P1 ELEMENT FOR LINEAR ELASTIC MATERIALS ---------------%
% ----------------------------------------------------------------------- %
% ------- By Dott. Ing. Marco Pingaro (P.h.D. Students)           ------- %                                   
% ------- University of Pavia                                     ------- %
% -------                                                         ------- %   
% ----------------------------------------------------------------------- %
%                                                                         %
% ------- Displacement  : P1                                      ------- %
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
poisson = input('Inserire modulo di Poisson = '); 
%-------------------------------------------------------------------------%
dx = length/ndx;
dy = height/ndy;
% Carico geometria dati 
[coordinates,nnod]=CoordinatesType2(ndx,ndy,dx,dy);
[element,nelem]=ElementType2(ndx,ndy);
[mc,ngdlu]=CorrispoMC(element,nelem,nnod);
%
%-------------------------------------------------------------------------%
% GLOBAL STIFF MATRIX
A = sparse(ngdlu,ngdlu);
b = sparse(ngdlu,1);
% ASSEMBLY GLOBAL MATRIX A, B and b (Global Load Vector)  
for k = 1:nelem
    P([1 2])=coordinates(element(k,1),[1 2]);
    P([3 4])=coordinates(element(k,2),[1 2]);
    P([5 6])=coordinates(element(k,3),[1 2]);
    [AELEM,JF] = GradGrad(P,poisson,young);
    [bf] = BodyLoad(JF,f);
    for i=1:6
        for j=1:6
           A(mc(k,i),mc(k,j)) = A(mc(k,i),mc(k,j))+AELEM(i,j);
        end
        b(mc(k,i),1) = b(mc(k,i),1) + bf(i,1);
    end
end
%-------------------------------------------------------------------------%
% Vincoli EX : 
%      Trave incernierata e incarrellata ai due estremi bassi
%-------------------------------------------------------------------------%
I = [1, 2, 2*(ndx+1)];
for jj = 1:size(I,2)
     A(I(jj),:) = 0;
     A(I(jj),I(jj)) = 1;
     b(I(jj),1) = 0;
end
%-------------------------------------------------------------------------%
% Solution
soluz = A\b;
Ux = soluz(1:2:ngdlu); 
Uy = soluz(2:2:ngdlu);
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