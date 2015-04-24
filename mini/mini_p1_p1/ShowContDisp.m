function ShowContDisp(element,coordinates,r)
hold on
for j=1:size(element,1)
    u(1,[1 2 3]) = r(element(j,[1 2 3]),1);
    trisurf([1 2 3],coordinates(element(j,:),1),...
        coordinates(element(j,:),2),...
        u,'facecolor','interp');
end
axis auto
view(-60,50);
