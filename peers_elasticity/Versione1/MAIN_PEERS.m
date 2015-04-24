%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        MAIN OF PEERS ELEMENT                            %       
%  by Dott. Ing. Marco Pingaro (PhD student)                              %
%  Department of Civil Engineering and Architecture,                      %
%             Via Ferrata 3, 27100, Pavia                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  - Weak Formulation of the Poblem :                                     %
%          < C*Sig,Tau > + < u, div(Tau) > + < u, as(Tau) >= < u, Tau*n > %
%                          < v, div(Sig) >                 = < f, v >     %
%                                            < v, as(Sig) >= 0            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
%% GEOMETRY AND PARAMETERS
lunghezza = 1; altezza = 1;
ndx = 10; ndy = 10;
dx = lunghezza/ndx; dy = altezza/ndy;
f(1,1) = 1; f(2,1) = 1;      % Carico per unita di superficie.
rho = 1;                     % Densita materiale
young = 1;                   % Costanti Ingegneristiche 
poisson = 0.3;               % Costanti Ingegneristiche.
tstart = tic;
%%
[coordinates]=CoordinatesType1(ndx,ndy,dx,dy);
[element]=ElementType1(ndx,ndy);
[mc]=CorrispoType1(ndx,ndy);
[nodes2element,nodes2edge,noedges,edge2element,...
exterioredge]=edge(element,coordinates);

nelem = size(element,1);     % Numero complessivo di elementi.
nnod = size(coordinates,1);  % Numero di nodi
ngdlu = 2*nelem;             % Numero gradi di liberta  spostamenti.
ngdls = max(max(mc));        % Numero gradi di liberta  sforzo.
ngdlse = 8;                  % Numero gradi di liberta  sforzo elementare. 
ngdltot = ngdls+ngdlu+nnod;  % Numero gradi di liberta  intero sistema.
lambda = young*poisson/((1+poisson)*(1-2*poisson));
mu = young/(2*(1+poisson));

%% GLOBAL STIFF MATRIX
A = sparse(ngdls,ngdls); 
B = sparse(ngdls,2*nelem);
C = sparse(ngdls,nnod);
M = sparse(ngdlu,ngdlu);
bf = sparse(2*nelem,1);

%% ASSEMBLY GLOBAL MATRIX A B & C
for k=1:nelem
    P([1 2],1)=coordinates(element(k,1),[1 2]);
    P([3 4],1)=coordinates(element(k,2),[1 2]);
    P([5 6],1)=coordinates(element(k,3),[1 2]);
    s = mc(k,ngdlse+1);
    [AELEM,BELEM,CELEM,AREA] = SigmaTauPEERS(P,mu,lambda,s);
    MELEM = mass_matrix(P,rho);
    bf(2*(k-1)+[1 2],1) = [f(1,1)*AREA, f(2,1)*AREA];
    for i=1:ngdlse
        for j=1:ngdlse
            if i>j
                AELEM(i,j)=AELEM(j,i);
            end
            A(mc(k,i),mc(k,j)) = A(mc(k,i),mc(k,j))+AELEM(i,j);
        end
        B(mc(k,i),2*(k-1)+[1 2]) = BELEM(i,[1 2]);
        for z =1:3
            C(mc(k,i), element(k,z)) = C(mc(k,i),element(k,z))...
                +CELEM(i,z);
        end
    end
    M(2*(k-1)+1, 2*(k-1)+1) = MELEM(1,1);
    M(2*(k-1)+2, 2*(k-1)+2) = MELEM(2,2);
end
K = [A, B, C; B', sparse(2*nelem,2*nelem+nnod); C', sparse(nnod,2*nelem+nnod)];
% GLOBAL VECTOR LOAD
bg = sparse(ngdls,1);
bh = sparse(nnod,1);
p = [bg;bf;bh]; 

clear i j k s lambda mu
clear AELEM BELEM CELEM 
clear A B C P AREA
clear bf bg bh
%----------------------------------------------------------------
%[p] = Load(f,ngdls,nelem);
%% IMPOSITION OF THE BOUNDARY CONDITION (Neumann)
% u_D == 0 non faccio niente
% V(1) == 0 or 1 incastro o svincolo lato orizzontale basso
% V(2) == 0 or 1 incastro o svoncolo lato verticale destro
% V(3) == 0 or 1 incastro o svincolo lato orizzontale alto
% V(4) == 0 or 1 incastro o svincolo lato verticale sinistro
%V =[0 0 0 0];
%[K,p] = NeumannBoundary(ndx,ndy,K,p,V);
flux1 = 1:ndx*2; 
flux2 = 6*ndx+1:1:3*nelem+2*ndy-1;
flux3 = 3*nelem+2*ndy+1:3*nelem+2*ndx+2*ndy;
flux4 = [];

%V= [flux1,flux3];
V = [];
nnv = size(V,2);
for i = 1:nnv
    K(V(i),:) = 0;
    K(V(i),V(i))= 1;
    p(V(i),1) = 0;
end

%% SOLUTION OF THE SYSTEM
soluz = K\p;

%% POST PROCESSING
stress = soluz(1:ngdls);
etasp = soluz(ngdls+ngdlu+1:ngdltot);
spost = soluz(ngdls+1:ngdls+ngdlu,1);
ux = spost(1:2:ngdlu); uy = spost(2:2:ngdlu);
defo = DefoMesh(coordinates, ux, uy, nodes2element);
ValStr = StressValue(mc,stress,nelem);
tfinish = toc(tstart);
%% PLOT
% Define structure solution and mesh
msh.coordinates = coordinates;
msh.element     = element;
solu.defo      = defo;
solu.ux        = ux;
solu.uy        = uy;
solu.etasp     = etasp;
solu.ValStr    = ValStr;

plt_solution(msh,solu)
%% END OF PROGRAM