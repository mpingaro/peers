function ShowRotRigid(element,coordinates,rot)
hold on
for j=1:size(element,1)
    u(1,[1 2 3]) = rot(element(j,[1 2 3]),1);
    trisurf([1 2 3],coordinates(element(j,:),1),...
        coordinates(element(j,:),2),...
        u,'facecolor','interp');
end 
title('Displacement Rotation')
view(-60,50);
