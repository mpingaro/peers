function plotmeshcook(EL,COORDINATES,NX,NY,TYPE)
if TYPE == 1
    figure, trimesh (EL,COORDINATES(:,1),COORDINATES(:,2))
    title 'UNDEFORMED MESH'
    xlabel 'Beam length'
    ylabel 'Beam height'
    axis equal
elseif TYPE == 2
    X = reshape(COORDINATES(:,1),NX+1,NY+1);
    Y = reshape(COORDINATES(:,2),NX+1,NY+1);
    figure, mesh(X,Y,zeros(NX+1,NY+1))
    title('UNDEFORMED MESH')
    xlabel 'Beam length'
    ylabel 'Beam height'
    view(0,90)
    axis equal
    grid off
end
end

