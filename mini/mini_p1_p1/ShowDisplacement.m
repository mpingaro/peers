function ShowDisplacement(element,coordinates,u)
hold on
for j=1:size(element,1)
    trisurf([1 2 3],coordinates(element(j,:),1),...
        coordinates(element(j,:),2),...
        ones(3,1)*u(j)','facecolor','interp');
end 
title('Displacement')
axis auto
view(-60,50);
