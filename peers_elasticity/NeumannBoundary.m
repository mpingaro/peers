function [flux1n,flux1t,flux2n,flux2t,flux3n,flux3t,flux4n,flux4t] =...
     NeumannBoundary(ndx,ndy)
% Numerazione flussi di lato per la mesh tipo 1:

% Flussi lato 1 (orizzontale basso)
flux1n = [1:ndx*2]; 
flux1t = [];
% Flussi lato 2 (verticale destro)
flux2n = [6*ndx+1:1:3*nel+2*ndy-1];
flux2t = [];
% Flussi lato 3 (orizzontale alto)
flux3n = [3*nel+2*ndy+1:3*nel+2*ndx+2*ndy];
flux3t = [];
% Flussi lato 4 (verticale sinistro)
flux4n = [];
flux4t = [];

return
