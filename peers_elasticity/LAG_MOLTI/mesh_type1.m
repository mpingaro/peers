function geo = mesh_type1(ndx,ndy,dx,dy)

% Coordinates
npoint = (ndx+1)*(ndy+1);
coordinates = zeros(npoint,1);
l = ndx*dx; h = ndy*dy;
vecx = 0:dx:l; vecy = 0:dy:h;
for i=1:ndx+1
    y(i,:) = vecy;
end
for j=1:ndy+1
   x(:,j) =vecx; 
end
coordinates(:,1) = reshape(x,npoint,1);
coordinates(:,2) = reshape(y,npoint,1);

% Elements
nelem = ndx*ndy*2;
element = zeros(nelem,3);

for i = 1:ndy
    for j = 1:ndx
        i1 = (ndx+1)*(i-1)+j;
        i2 = i1+1;
        i3 = (ndx+1)*i+j;
        i4 = i3+1;
        
        a = 2*ndx*(i-1)+2*(j-1)+1;
        b = a+1;
        element(a,[1 2 3]) = [i1, i2, i3];
        element(b,[1 2 3]) = [i4, i3, i2];
    end
end

% Matrix nodes2element
nodes2element=sparse(npoint,npoint);
for j=1:nelem
    nodes2element(element(j,:),element(j,[2 3 1]))= ...
        nodes2element(element(j,:),element(j,[2 3 1]))+j*eye(3,3);
end
% Matrix nodes2edge
B=nodes2element+nodes2element'; [I,J]=find(triu(B));
nodes2edge=sparse(I,J,1:size(I,1),npoint,npoint);
nodes2edge=nodes2edge+nodes2edge';

% Noedges 
noedges=size(I,1);

%Matrix edge2element
edge2element=zeros(noedges,4);
for m=1:nelem
    for k=1:3
        initial_node = element(m,k);
        end_node=element(m,rem(k,3)+1);
        p=nodes2edge(element(m,k),element(m,rem(k,3)+1));
        if edge2element(p,1)==0
            edge2element(p,:)=[initial_node end_node ...
                nodes2element(initial_node,end_node) ...
                nodes2element(end_node,initial_node)];
        end
    end
end

interioredge=find(edge2element(:,4));
exterioredge=find(edge2element(:,4)==0);

% Struttura geo
geo.coordinates = coordinates;
geo.element = element;
geo.nod2elm = nodes2element;
geo.nod2edg = nodes2edge;
geo.edg2elm = edge2element;
geo.noedges = noedges;
geo.intedge = interioredge;
geo.extedge = exterioredge;
geo.nnod = npoint;
geo.nelem = nelem;


return