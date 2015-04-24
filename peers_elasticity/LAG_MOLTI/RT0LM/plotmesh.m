function plotmesh(geo)

title('Triangular Mesh','fontsize',14);
trimesh(geo.element,geo.coordinates(:,1),geo.coordinates(:,2),'color','g')
axis equal

end