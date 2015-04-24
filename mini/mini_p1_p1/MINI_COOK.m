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
% Carico geometria dati
sizemesh = input('Inserire la dimensione della mesh (Valore compreso tra 0 e 16) = ');
f(1,1) = input('Inserire forzante di superficie direzione x = '); 
f(2,1) = input('Inserire forzante di superficie direzione y = ');
g(1,1) = input('Inserire forzante sul lato destro direzione (distribuito) x = ');
g(2,1) = input('Inserire forzante sul lato destro direzione (distribuito) y = ');
young  = input('Inserire modulo di Young = ');

%-------------------------------------------------------------------------%
pv =[0.0, 0.0; 48.0, 44.0; 48.0, 60.0; 0.0, 44.0; 0.0, 0.0];
[coordinates,element] = ...
    distmesh2d(@dpoly,@huniform,sizemesh,[0,0; 48,60],pv,pv);
[nodes2element,nodes2edge,noedges,edge2element,...
exterioredge]=edge(element,coordinates);
%-------------------------------------------------------------------------%
nnod = size(coordinates,1);
nelem = size(element,1);
[mc,ngdlu]=CorrispoMC(element,nelem,nnod);
%
ngdlp = nnod;
ngdlt = ngdlu+ngdlp;
mu2 = young/1.5;
%-------------------------------------------------------------------------%
% GLOBAL STIFF MATRIX
A =  sparse(ngdlu,ngdlu);
B =  sparse(ngdlu,ngdlp);
bf = sparse(ngdlu,1);
% ASSEMBLY GLOBAL MATRIX A, B and b (Global Load Vector)  
for k = 1:nelem
    P([1 2])=coordinates(element(k,1),[1 2]);
    P([3 4])=coordinates(element(k,2),[1 2]);
    P([5 6])=coordinates(element(k,3),[1 2]);
    [AELEM,BELEM,JF] = GradGrad(P);
    [lf] = BodyLoad(JF,f);
    for i=1:8
        for j=1:8
           A(mc(k,i),mc(k,j)) = A(mc(k,i),mc(k,j))+mu2*AELEM(i,j);
        end
        bf(mc(k,i),1) = bf(mc(k,i),1) + lf(i,1);
        for ii = 1:3
            B(mc(k,i),element(k,ii)) = B(mc(k,i),element(k,ii))...
                +BELEM(i,ii);
        end
    end
end
K = [A, -B; B', sparse(nnod,nnod)];
clear A B AELEM BELEM JF P
clear young poisson f
%-------------------------------------------------------------------------%
% Vincoli EX : 
%      Mensola di Cook's incastrata lato sinistro
%-------------------------------------------------------------------------%
[UD,UF] = BDcooks(coordinates,edge2element,exterioredge);
[bg] = BNcooks(coordinates,UF,g,ngdlu);
%
for jj = 1:size(UD,1)
     I = 2*UD(jj,1);
     II= 2*UD(jj,1)-1;
     K(I,:)  = 0;
     K(II,:) = 0;
     K(I,I)  = 1;
     K(II,II)= 1;
     bf(I,1) = 0;
     bf(II,1)= 0;
end
bt = bf + bg;
b = [bt; sparse(ngdlp,1)]; 
% ------------------------------------------------------------------------%
% SOLUTION
%
soluz = K\b;
Ux = soluz(1:2:2*nnod); 
Uy = soluz(2:2:2*nnod);
pressure = soluz(ngdlu+1:ngdlt);
defo =coordinates + [Ux, Uy];
%
%maxUy = max(Uy);
%fprintf('Massimo spostamanto verticale Point B = %d \n', maxUy);
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