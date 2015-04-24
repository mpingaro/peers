%==========================================================================
%                  C O O K ' S    C A N T I L E V E R
%                 Two-dimensional mixed finite element
%                           (by Marco Pingaro)
%==========================================================================
clear all; 
close all; 
clc

%% INPUT DATA
[NX,NY,TYPE] = inputcook();
NODES = [0 0; 48 44; 48 60; 0 44]; 
%NODES = [0 0; 50 20; 50 40; 0 60];
DL1 = NODES(3,2)-NODES(2,2); 
DL2 = NODES(4,2);

%% COORDINATES MATRIX
[COORDINATES] = coordcook(NODES,NX,NY,DL1,DL2);

%% ELEMENT MATRIX
[EL] = elcook(NX,NY,TYPE);

%% CORRESPONDENCE MATRIX (element)
[MCORR] = mccook(EL,TYPE);

%% CORRESPONDENCE MATRIX (stress flux)
[FLUX] = fluxcook(NX,NY,TYPE);

%% BOUNDARY CONDITIONS
BC = [];
for i = 0:NY
    BC = [BC , 1+i*(NX+1)];    
end

%% PLOT
plotmeshcook(EL,COORDINATES,NX,NY,TYPE)