function [mc,ngdlu]=CorrispoMC(element,nelem,nnod,ndx,ndy)

mc = zeros(nelem,12);         % prealloco matrice di corrispondenza
ngdlu = 2*nnod;

mc(:,1) = 2*element(:,1)-1;
mc(:,2) = 2*element(:,1);
mc(:,3) = 2*element(:,2)-1;
mc(:,4) = 2*element(:,2);
mc(:,5) = 2*element(:,3)-1;
mc(:,6) = 2*element(:,3);

for i = 1:ndy
    for j = 1:ndx

        i1 = ngdlu+(6*ndx+2)*(i-1)+2*(j-1)+1;
        i2 = ngdlu+(6*ndx+2)*(i-1)+2*ndx+4*(j-1)+1;
        i3 = i2+4;
        i4 = ngdlu+6*ndx+2+(6*ndx+2)*(i-1)+2*(j-1)+1;
        i5 = i2+2;
        
        a = 2*ndx*(i-1) + 2*(j-1)+1;
        b = a+1;
        
        mc(a,7) = i1;
        mc(a,8) = i1+1;
        mc(a,9) = i5;
        mc(a,10)= i5+1;
        mc(a,11)= i2;
        mc(a,12)= i2+1;
        
        mc(b,7) = i4;
        mc(b,8) = i4+1;
        mc(b,9) = i5;
        mc(b,10)= i5+1;
        mc(b,11)= i3;
        mc(b,12)= i3+1;
    end
end


return