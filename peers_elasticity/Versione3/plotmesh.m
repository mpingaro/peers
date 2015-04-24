function figura = plotmesh(element,coordinates)
axis equal 
%hold on
title('Triangular Mesh','fontsize',14);
trimesh(element,coordinates(:,1),coordinates(:,2),'color','g')

return