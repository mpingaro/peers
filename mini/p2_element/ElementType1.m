function [element,nelem]=ElementType1(ndx,ndy)

nelem = ndx*ndy*2;
element = zeros(nelem,3);

for i = 1:ndy
    for j = 1:ndx
        i1 = (ndx+1)*(i-1)+j;
        i2 = i1+1;
        i3 = (ndx+1)*i+j;
        i4 = i3+1;
        
        a = 2*ndx*(i-1)+2*(j-1)+1;
        b = a+1;
        element(a,[1 2 3]) = [i1, i2, i3];
        element(b,[1 2 3]) = [i4, i3, i2];
    end
end

return