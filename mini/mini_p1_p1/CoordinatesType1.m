function [coordinates,nnod]=CoordinatesType1(ndx,ndy,dx,dy)

nnod = (ndx+1)*(ndy+1);
coordinates = zeros(nnod,1);
l = ndx*dx; h = ndy*dy;
vecx = [0:dx:l]; vecy = [0:dy:h];
for i=1:ndx+1
    y(i,:) = vecy;
end
for j=1:ndy+1
   x(:,j) =vecx; 
end
coordinates(:,1) = reshape(x,nnod,1);
coordinates(:,2) = reshape(y,nnod,1);

return