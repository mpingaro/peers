function [nodes2element,nodes2edge,noedges,edge2element,...
exterioredge] = edge(geo)

% Matrix nodes2element
nodes2element=sparse(size(geo.coordinates,1),size(geo.coordinates,1));
for j=1:size(geo.element,1)
    nodes2element(geo.element(j,:),geo.element(j,[2 3 1]))= ...
        nodes2element(geo.element(j,:),geo.element(j,[2 3 1]))+j*eye(3,3);
end
% Matrix nodes2edge
B=nodes2element+nodes2element'; [I,J]=find(triu(B));
nodes2edge=sparse(I,J,1:size(I,1),size(geo.coordinates,1),size(geo.coordinates,1));
nodes2edge=nodes2edge+nodes2edge';

% Noedges 
noedges=size(I,1);

%Matrix edge2element
edge2element=zeros(noedges,4);
for m=1:size(geo.element,1)
    for k=1:3
        initial_node=geo.element(m,k);
        end_node=geo.element(m,rem(k,3)+1);
        p=nodes2edge(geo.element(m,k),geo.element(m,rem(k,3)+1));
        if edge2element(p,1)==0
            edge2element(p,:)=[initial_node end_node ...
                nodes2element(initial_node,end_node) ...
                nodes2element(end_node,initial_node)];
        end
    end
end

%interioredge=find(edge2element(:,4));
exterioredge=find(edge2element(:,4)==0);

return

