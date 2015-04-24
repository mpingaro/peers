function [coordinates]=CoordinatesType1(ndx,ndy,dx,dy)

npoint = (ndx+1)*(ndy+1);
coordinates = zeros(npoint,1);
l = ndx*dx; h = ndy*dy;
vecx = [0:dx:l]; vecy = [0:dy:h];
for i=1:ndx+1
    y(i,:) = vecy;
end
for j=1:ndy+1
   x(:,j) =vecx; 
end
coordinates(:,1) = reshape(x,npoint,1);
coordinates(:,2) = reshape(y,npoint,1);
%save coordinates
return