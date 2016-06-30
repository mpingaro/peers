%-------------------------------------------------------------------------%
%----------------   PEERS ELEMENT IN LINEAR ELASTICITY -------------------%
%----------------                                      -------------------% 
%  Programming By Dott. Ing. Marco Pingaro (University of Pavia)          %
%-------------------------------------------------------------------------%
clear all; close all; clc
%-- GEOMETRY OF BEAM   ---------------------------------------------------%
lx = input('Length in direction x of the beam = '); 
ly = input('Length in direction y of the beam = ');
% Number Of Suddivision
ndx = input('Number of suddivision in direction x = '); 
ndy = input('Number of suddivision in direction y = ');
% Elastic Modulus
young = 1;                                             % Young Modulus 
poisson = 0.3;                                         % Poisson Modulus
% Body Forces (for unit surface)
f(1,1) = input('Body forces in direction x = ');
f(2,1) = input('Body forces in direction y = ');
tstart = tic;
%-------------------------------------------------------------------------%
dx = lx/ndx; dy = ly/ndy;
% Prepare Input Geometry Data
[coordinates]=CoordinatesType1(ndx,ndy,dx,dy);
[element]=ElementType1(ndx,ndy);
[mc]=CorrispoType1(ndx,ndy);
%-------------------------------------------------------------------------%
% Vector Of Principal Index
ngdls = max(max(mc));
nelem = size(element,1);
nnod = size(coordinates,1);
ngdlu = 2*nelem;
ngdltot = ngdls+ngdlu+nnod;

lambda = young*poisson/((1+poisson)*(1-2*poisson));
mu = young/(2*(1+poisson));

% GLOBAL STIFF MATRIX & GLOBAL VECTOR LOAD
[K,p] = assembly(coordinates,element,mc,mu,lambda,f);
%----------------------------------------------------------------
% IMPOSITION OF THE BOUNDARY CONDITION (Neumann)

% u_D == 0 non faccio niente
% V(1) == 0 or 1 incastro o svincolo lato orizzontale basso
% V(2) == 0 or 1 incastro o svoncolo lato verticale destro
% V(3) == 0 or 1 incastro o svincolo lato orizzontale alto
% V(4) == 0 or 1 incastro o svincolo lato verticale sinistro
%V =[0 0 0 0];
%[K,p] = NeumannBoundary(ndx,ndy,K,p,V);
flux1 = [1:ndx*2]; 
flux2 = [6*ndx+1:1:3*nelem+2*ndy-1];
flux3 = [3*nelem+2*ndy+1:3*nelem+2*ndx+2*ndy];
flux4 = [];
V= [flux1, flux3];
nnv = size(V,2);
for i = 1:nnv
    K(V(i),:) = 0;
    K(V(i),V(i))= 1;
    p(V(i),1) = 0;
end
%----------------------------------------------------------------
% SOLUTION OF THE SYSTEM
soluz = K\p;
etasp = soluz(ngdls+ngdlu+1:ngdltot);
spost = soluz(ngdls+1:ngdls+ngdlu,1);
ux = spost(1:2:ngdlu); uy = spost(2:2:ngdlu);

tfinish = toc(tstart)
% PLOT
figure
plotmesh(element,coordinates)
figure
ShowDisplacement(element,coordinates,ux)
title('Displacement Ux')
figure
ShowDisplacement(element,coordinates,uy)
title('Displacement Uy')
figure
ShowRotRigid(element,coordinates,etasp)
