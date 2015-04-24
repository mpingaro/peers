function [EELEM] = PEERS_intcontinuity(geo,nedge)

EELEM = zeros(8,2,2);

I = geo.edg2elm(nedge,:);

% Baricientri triangoli adiacenti lati
bari1 = sum( geo.coordinates(geo.element(I(1,3),[1 2 3]),:) )/3;  
bari2 = sum( geo.coordinates(geo.element(I(1,4),[1 2 3]),:) )/3;
% Lato
EV = geo.coordinates(I(1,2),:)-geo.coordinates(I(1,1),:);
gamma1 = EV(2)*( geo.coordinates(I(1,1),1)-bari1(1) ) - EV(1)*( geo.coordinates(I(1,1),2)-bari1(2) );
gamma2 = EV(2)*( geo.coordinates(I(1,1),1)-bari2(1) ) - EV(1)*( geo.coordinates(I(1,1),2)-bari2(2) );
% 
EELEM([1 2 3 4 5 6],1,1) = [-EV(2), -EV(2), EV(1), -EV(1), -gamma1, -gamma1]; 
EELEM([1 2 3 4 5 6],2,1) = [EV(1), -EV(1), -EV(2), -EV(2), -gamma1,  gamma1];

EELEM([1 2 3 4 5 6],1,2) = [EV(2), EV(2), -EV(1), EV(1), gamma2,  gamma2]; 
EELEM([1 2 3 4 5 6],2,2) = [-EV(1), EV(1), EV(2), EV(2), gamma2, -gamma2];
end