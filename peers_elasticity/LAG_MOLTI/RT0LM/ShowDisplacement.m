function ShowDisplacement(geo,u)
hold on
for j=1:geo.nelem
    trisurf([1 2 3],geo.coordinates(geo.element(j,:),1),...
        geo.coordinates(geo.element(j,:),2),...
        ones(3,1)*u.disp(j)','facecolor','interp');
end 
title('Displacement')
axis auto
view(-60,50);

end