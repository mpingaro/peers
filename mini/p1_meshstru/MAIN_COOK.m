% ----------------------------------------------------------------------- %
% ---------------- P1 ELEMENT FOR LINEAR ELASTIC MATERIALS ---------------%
% ----------------------------------------------------------------------- %
% ------- By Dott. Ing. Marco Pingaro (P.h.D. Student)            ------- %                                   
% ------- University of Pavia                                     ------- %
% -------                                                         ------- %   
% ----------------------------------------------------------------------- %
%                                                                         %
% ------- Displacement  : P1 or Q1 element                        ------- %
% ------- Cook's membrane problem                                 ------- %
% ----------------------------------------------------------------------- %
clear all; close all; clc;
%% INPUT
%
[NX,NY,TYPE] = inputcook();
NODES = [0 0; 48 44; 48 60; 0 44]; 
DL1 = NODES(3,2)-NODES(2,2); 
DL2 = NODES(4,2);
% COORDINATES MATRIX
[coordinates] = coordcook(NODES,NX,NY,DL1,DL2);
% ELEMENT MATRIX
[element] = elcook(NX,NY,TYPE);
% CORRISPONDENCE MATRIX
[mc] = mccook(element,TYPE);
%
f(1,1)  = input('Inserire forzante di superficie direzione x = '); 
f(2,1)  = input('Inserire forzante di superficie direzione y = ');
f_traction  = input('Inserire forzante sul lato destro direzione y = ');
young   = input('Inserire modulo di Young = ');
poisson = input('Inserire modulo di Poisson = '); 
%-------------------------------------------------------------------------%
nnod = size(coordinates,1); 
nelem = size(element,1);
ngdlu = 2*nnod;
ngdlel = size(mc,2);
%-------------------------------------------------------------------------%
% GLOBAL STIFF MATRIX
A = sparse(ngdlu,ngdlu);
bf = sparse(ngdlu,1);
% ASSEMBLY GLOBAL MATRIX A, B and b (Global Load Vector)  
for k = 1:nelem
    if TYPE == 1
        P([1 2],1)=coordinates(element(k,1),[1 2]);
        P([3 4],1)=coordinates(element(k,2),[1 2]);
        P([5 6],1)=coordinates(element(k,3),[1 2]);
        [AELEM,JF] = GradGrad(P,poisson,young);
        [lf] = BodyLoad(JF,f);
    elseif TYPE == 2
        P(1,[1 2])=coordinates(element(k,1),[1 2]);
        P(2,[1 2])=coordinates(element(k,2),[1 2]);
        P(3,[1 2])=coordinates(element(k,3),[1 2]);
        P(4,[1 2])=coordinates(element(k,4),[1 2]);
        [AELEM,lf] = GradGrad_Quad(P,young,poisson,f);
    end
    for i=1:ngdlel
        for j=1:ngdlel
           A(mc(k,i),mc(k,j)) = A(mc(k,i),mc(k,j))+AELEM(i,j);
        end
        bf(mc(k,i),1) = bf(mc(k,i),1) + lf(i,1);
    end
end
%-------------------------------------------------------------------------%
% Vincoli EX : 
%      Mensola di Cook's incastrata lato sinistro
%-------------------------------------------------------------------------%
%% BOUNDARY IMPOSITION 
%     Trave incastrata e trazione sul lato verticale destro

% Calcoli indici bordo incastrato e caricato
tra = zeros(NY+1,1);
inc = zeros(2*(NY+1),1);
for i= 1:NY+1
    tra(i,1) = 2*(NX+1)*i;
    inc(2*(i-1)+1,1) = 2*((NX+1)*(i-1)+1)-1;
    inc(2*(i-1)+2,1) = 2*((NX+1)*(i-1)+1);
end
% Imposizione vincoli
for I = 1:size(inc,1)
    A(inc(I,1), :) = 0;
    A(inc(I,1),inc(I,1)) = 1;
    bf(inc(I,1),1) = 0;
end

f_traction = f_traction/(2*NY);
for J = 1:size(tra,1)
    if J==1 || J == size(tra,1)
        bg(tra(J,1),1) = f_traction;
    else
        bg(tra(J,1),1) = 2*f_traction;
    end
end
b = bf + bg;
%-------------------------------------------------------------------------%
% SOLUTION
soluz = A\b;
Ux = soluz(1:2:ngdlu-1); 
Uy = soluz(2:2:ngdlu);
defo =coordinates + [Ux, Uy];

% PLOT SOLUTION
plotmeshcook(element,coordinates,NX,NY,TYPE)
title('Undeformed Mesh','fontsize',14);
%
plotmeshcook(element,defo,NX,NY,TYPE)
title('Deformed Mesh','fontsize',14);