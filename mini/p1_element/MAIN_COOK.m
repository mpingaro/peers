% ----------------------------------------------------------------------- %
% ---------------- P1 ELEMENT FOR LINEAR ELASTIC MATERIALS ---------------%
% ----------------------------------------------------------------------- %
% ------- By Dott. Ing. Marco Pingaro (P.h.D. Student)            ------- %                                   
% ------- University of Pavia                                     ------- %
% -------                                                         ------- %   
% ----------------------------------------------------------------------- %
%                                                                         %
% ------- Displacement  : P1                                      ------- %
% ------- Cook's membrane problem                                 ------- %
% ----------------------------------------------------------------------- %
clear all; close all; clc;
%% INPUT
% INPUT DATA
[NX,NY,TYPE] = inputcook();
NODES = [0 0; 48 44; 48 60; 0 44]; 
DL1 = NODES(3,2)-NODES(2,2); 
DL2 = NODES(4,2);
% COORDINATES MATRIX
[coordinates] = coordcook(NODES,NX,NY,DL1,DL2);
% ELEMENT MATRIX
[element] = elcook(NX,NY,TYPE);
%
[nodes2element,nodes2edge,noedges,edge2element,...
exterioredge]=edge(element,coordinates);
f(1,1)  = input('Inserire forzante di superficie direzione x = '); 
f(2,1)  = input('Inserire forzante di superficie direzione y = ');
g(1,1)  = input('Inserire forzante sul lato destro direzione x = ');
g(2,1)  = input('Inserire forzante sul lato destro direzione y = ');
young   = input('Inserire modulo di Young = ');
poisson = input('Inserire modulo di Poisson = '); 
%-------------------------------------------------------------------------%
nnod = size(coordinates,1); nelem = size(element,1);
[mc,ngdlu]=CorrispoMC(element,nelem,nnod);
%-------------------------------------------------------------------------%
% GLOBAL STIFF MATRIX
A = sparse(ngdlu,ngdlu);
bf = sparse(ngdlu,1);
% ASSEMBLY GLOBAL MATRIX A, B and b (Global Load Vector)  
for k = 1:nelem
    P([1 2],1)=coordinates(element(k,1),[1 2]);
    P([3 4],1)=coordinates(element(k,2),[1 2]);
    P([5 6],1)=coordinates(element(k,3),[1 2]);
    [AELEM,JF] = GradGrad(P,poisson,young);
    [lf] = BodyLoad(JF,f);
    for i=1:6
        for j=1:6
           A(mc(k,i),mc(k,j)) = A(mc(k,i),mc(k,j))+AELEM(i,j);
        end
        bf(mc(k,i),1) = bf(mc(k,i),1) + lf(i,1);
    end
end
%-------------------------------------------------------------------------%
% Vincoli EX : 
%      Mensola di Cook's incastrata lato sinistro
%-------------------------------------------------------------------------%
[UD,UF] = BDcooks(coordinates,edge2element,exterioredge);
[bg] = BNcooks(coordinates,UF,g,ngdlu);
% Imposizione delle condizioni di Dirichlet
for jj = 1:size(UD,1)
     I = 2*UD(jj,1);
     II= 2*UD(jj,1)-1;
     A(I,:)  = 0;
     A(II,:) = 0;
     A(I,I)  = 1;
     A(II,II)= 1;
     bf(I,1) = 0;
     bf(II,1)= 0;
end
b = bf + bg;
%-------------------------------------------------------------------------%
% SOLUTION
soluz = A\b;
Ux = soluz(1:2:ngdlu); 
Uy = soluz(2:2:ngdlu);
defo =coordinates + [Ux, Uy];

% PLOT SOLUTION
figure
plotmesh(element,coordinates)
title('Undeformed Mesh','fontsize',14);
figure
plotmesh(element,defo)
title('Deformed Mesh','fontsize',14);