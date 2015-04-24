function [bg] = BNcooks(coordinates,UF,g,ngdlu)

bg = sparse(ngdlu,1);

for i = 1:size(UF,1)
    lato = sqrt((coordinates(UF(i,1),1)-coordinates(UF(i,2),1))^2 +...
        (coordinates(UF(i,1),2)-coordinates(UF(i,2),2))^2);
    gnodex = g(1,1)*lato/2;
    gnodey = g(2,1)*lato/2;
    bg(2*UF(i,1)-1,1) = bg(2*UF(i,1)-1,1) + gnodex;
    bg(2*UF(i,1),1)   = bg(2*UF(i,1),1) + gnodey;
    bg(2*UF(i,2)-1,1) = bg(2*UF(i,2)-1,1) + gnodex;
    bg(2*UF(i,2),1)   = bg(2*UF(i,2),1) + gnodey;
end

return