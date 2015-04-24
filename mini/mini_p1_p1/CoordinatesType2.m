function [coordinates,nnod]=CoordinatesType2(ndx,ndy,dx,dy)

n1 = (ndx+1)*(ndy+1);
n2 = ndx*ndy;
nnod = n1+n2;
coordinates = zeros(nnod,1);
l = ndx*dx; h = ndy*dy;
vecx = [0:dx:l]; vecy = [0:dy:h];
intx = [dx/2:dx:l]; inty = [dy/2:dy:h];

for i=1:ndx+1 y(i,:) = vecy; end
for i=1:ndx yi(i,:) = inty; end
for j=1:ndy+1 x(:,j) =vecx; end
for j=1:ndy xi(:,j) = intx; end

cor1(:,1) = reshape(x,n1,1);
cor1(:,2) = reshape(y,n1,1);
cor2(:,1) = reshape(xi,n2,1);
cor2(:,2) = reshape(yi,n2,1);

coordinates = [cor1;cor2];

return