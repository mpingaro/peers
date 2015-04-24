function ShowStress(element,coordinates,str)
hold on
for j=1:size(element,1)
    u(1,[1 2 3]) = str(3*(j-1)+[1 2 3], 1);
    trisurf([1 2 3],coordinates(element(j,:),1),...
        coordinates(element(j,:),2),...
        u,'facecolor','interp');
end 
view(0,90);
axis equal