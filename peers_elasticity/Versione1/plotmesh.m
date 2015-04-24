function figura = plotmesh(element,coordinates)

title('Triangular Mesh','fontsize',14);
trimesh(element,coordinates(:,1),coordinates(:,2),'color','g')
axis equal
return