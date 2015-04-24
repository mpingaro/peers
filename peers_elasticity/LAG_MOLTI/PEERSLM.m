%----------------- PEERS VERSION LAGRANGIAN MULTIPLIER -------------------%
% -- By Dott. Ing. Marco Pingaro (PhD Student : IUSS Pavia, DICAr)        %
%                                                                         %
%-------------------------------------------------------------------------%
clear all; clc; close all;
% Elastic constant
young = 1;                   % Costanti Ingegneristiche 
poisson = 0.3;               % Costanti Ingegneristiche.
bf = [1,1];                  % Carico di superficie
%
ndx = 10;                     
ndy = 10; 
dx = 0.1; 
dy = 0.1;
%
lambda = young*poisson/((1+poisson)*(1-2*poisson)); % Costanti di Lamé
mu = young/(2*(1+poisson));                         % Costanti di Lamé
alpha = 1/(2*mu); 
beta = lambda/(4*mu*(mu+lambda));
% ----------------------------------------------------------------------- %
% Struttura geometria
geo = mesh_type1(ndx,ndy,dx,dy);
% Solve 
u = Solve_PEERS(geo,alpha,beta,bf);
% Geometria
figure, plotmesh(geo);
% Displacement Solution
figure, ShowDisplacement(geo,u.disp(1,:));
figure, ShowDisplacement(geo,u.disp(2,:));
% ----------------------------------------------------------------------- %
format long
disp( max(u.disp(1,:)) );
disp( min(u.disp(1,:)) );
%
disp(-0.003708369031603);
disp(-0.088201676737408);