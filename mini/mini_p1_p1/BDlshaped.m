function [UD,UF] = BDlshaped(coordinates,edge2element,exterioredge)

J=1;
JJ =1;
for i = 1:size(exterioredge,1)
     if (coordinates( edge2element(exterioredge(i,1),1),2)== 0)
         UD(J,1) = edge2element(exterioredge(i,1),1);
         J = J +1;
     end
     if (coordinates(edge2element(exterioredge(i,1),2),2) == 0)
         UD(J,1) = edge2element(exterioredge(i,1),2);
         J = J+1;
     end
     
     if (coordinates( edge2element(exterioredge(i,1),1),1) == 0 && ...
         coordinates( edge2element(exterioredge(i,1),2),1) == 0)
         UF(JJ,1) = edge2element(exterioredge(i,1),1);
         UF(JJ,2) = edge2element(exterioredge(i,1),2);
         JJ = JJ +1;
     end     
end


return