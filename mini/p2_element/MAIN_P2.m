%% --------------------------------------------------------------------- %%
% ---------------- P2 ELEMENT FOR LINEAR ELASTIC MATERIALS ---------------%
% ----------------------------------------------------------------------- %
% ------- By Dott. Ing. Marco Pingaro (P.h.D. Student)            ------- %                                   
% ------- University of Pavia                                     ------- %
% -------                                                         ------- %   
%                                                  xxx xxx        ------- %
%    MMM    MM    PPPPP                           x   x   x       ------- %
%    MMM    MM    PP   P                           xxx xxx        ------- %
%    MM MMMM M    PP   P                           V     V        ------- %
%    MM  MM  M    PPPPP                             V   V         ------- %
%    MM      M    PP        ------------------>      V V          ------- %         
%    MM      M o  PP    o   ------------------>       V           ------- %          
% ----------------------------------------------------------------------- %
%                                                                         %
% ------- Displacement  : P2                                      ------- %
% ----------------------------------------------------------------------- %
%% --------------------------------------------------------------------- %%
clear all; close all; clc;
%% -- INPUT ------------------------------------------------------------ %%

% --- PROVA
length = 4;
height = 2;
ndx = 4;
ndy = 2;
f(1,1) = 0;
f(2,1) = -1;
young = 70;
poisson = 0.3;

%length = input('Lunghezza della trave = ');
%height = input('Altezza della trave = ');
%ndx = input('Numero suddivisioni lungo x = ');
%ndy = input('Numero suddivisioni lungo y = ');
%f(1,1) = input('Inserire forzante di superficie direzione x = '); 
%f(2,1) = input('Inserire forzante di superficie direzione y = ');
%young = input('Inserire modulo di Young = ');
%poisson = input('Inserire modulo di Poisson = '); 
%% --------------------------------------------------------------------- %%
dx = length/ndx;
dy = height/ndy;
% Carico geometria dati 
[coordinates,nnod]=CoordinatesType1(ndx,ndy,dx,dy);
[element,nelem]=ElementType1(ndx,ndy);
[mc,ngdlu]=CorrispoMC(element,nelem,nnod,ndx,ndy);
ngdluu = 6*ndx*ndy+2*ndx+2*ndy;
ngdlt = ngdlu + ngdluu;
%
%% --------------------------------------------------------------------- %%
% GLOBAL STIFF MATRIX
A = sparse(ngdlt,ngdlt);
b = sparse(ngdlt,1);
% ASSEMBLY GLOBAL MATRIX A, B and b (Global Load Vector)  
for k = 1:nelem
    P(1,[1 2]) = coordinates(element(k,1),[1 2]);
    P(1,[3 4]) = coordinates(element(k,2),[1 2]);
    P(1,[5 6]) = coordinates(element(k,3),[1 2]);
    P(1,[7 8]) = [(P(1,1)+P(1,3))/2,(P(1,2)+P(1,4))/2];
    P(1,[9 10]) = [(P(1,3)+P(1,5))/2,(P(1,4)+P(1,6))/2];
    P(1,[11 12]) = [(P(1,5)+P(1,1))/2,(P(1,6)+P(1,2))/2];
    AELEM = p2element(P,poisson,young);
    bf = Load(f,P);
    for i=1:12
        for j=1:12
           A(mc(k,i),mc(k,j)) = A(mc(k,i),mc(k,j))+AELEM(i,j);
        end
        b(mc(k,i),1) = b(mc(k,i),1) + bf(i,1);
    end
end
%% --------------------------------------------------------------------- %%
% Vincoli EX : 
%      Trave incernierata e incarrellata ai due estremi bassi
% ----------------------------------------------------------------------- %
I = [1, 2, 2*(ndx+1), 2*(ndx+1)-1];
for jj = 1:size(I,2)
     A(I(jj),:) = 0;
     A(I(jj),I(jj)) = 1;
     b(I(jj),1) = 0;
end
%% --------------------------------------------------------------------- %%
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
axis equal