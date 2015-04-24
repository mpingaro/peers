function [DELEM] = RT0_intcontinuity(geo,nedge)

DELEM = zeros(3,1,2);

I = geo.edg2elm(nedge,:);

% Baricientri triangoli adiacenti lati
bari1 = sum( geo.coordinates(geo.element(I(1,3),[1 2 3]),:) )/3;  
bari2 = sum( geo.coordinates(geo.element(I(1,4),[1 2 3]),:) )/3;
% Lato
EV = geo.coordinates(I(1,2),:)-geo.coordinates(I(1,1),:);
gamma1 = EV(1)*( bari1(2) - geo.coordinates(I(1,1),2) ) + EV(2)*( geo.coordinates(I(1,1),1) - bari1(1) );
gamma2 = EV(1)*( bari2(2) - geo.coordinates(I(1,1),2) ) + EV(2)*( geo.coordinates(I(1,1),1) - bari2(1) );
% 
DELEM([1 2 3],1,1) = [-EV(2), EV(1), -gamma1]; 
DELEM([1 2 3],1,2) = [EV(2), -EV(1), gamma2]; 

end