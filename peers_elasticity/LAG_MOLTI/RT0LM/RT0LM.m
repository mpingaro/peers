%-------------- RT0 VERSION LAGRANGIAN MULTIPLIER ------------------------%
%----- by Dott. Ing. Marco Pingaro (PhD Student IUSS Pavia, DICAr) -------%
%-------------------------------------------------------------------------%
clear all; clc; close all;

bf = 1;                      % Carico di superficie
ndx = 2;                    % Suddivisioni in x                   
ndy = 2;                    % Suddivisioni in y
dx = 0.5;                   % Lunghezza elemento dir. x
dy = 0.5;                   % Lunghezza elemento dir. y

% Struttura geometria
geo = mesh_type1(ndx,ndy,dx,dy);
% Solve 
u = Solve_RT0(geo,bf);
% PLOT SOLUTION
% Geometria
figure, plotmesh(geo);
% Displacement Solution
figure, ShowDisplacement(geo,u);

format long
disp( max(u.disp(:,1)) );
disp( min(u.disp(:,1)) );