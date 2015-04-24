function [mc,ngdlu]=CorrispoMC(element,nelem,nnod)

mc = zeros(nelem,8);         % prealloco matrice di corrispondenza
ngdlu = 2*(nelem+nnod);

mc(:,1) = 2*element(:,1)-1;
mc(:,2) = 2*element(:,1);
mc(:,3) = 2*element(:,2)-1;
mc(:,4) = 2*element(:,2);
mc(:,5) = 2*element(:,3)-1;
mc(:,6) = 2*element(:,3);
mc(:,7) = [2*nnod+1:2:ngdlu];
mc(:,8) = [2*nnod+2:2:ngdlu];

return