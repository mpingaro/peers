function [element,nelem]=ElementType2(ndx,ndy)

nelem = ndx*ndy*4;
element = zeros(nelem,3);

for i = 1:ndy
    for j = 1:ndx
        i1 = (ndx+1)*(i-1)+j;
        i2 = i1+1;
        i3 = (ndx+1)*i+j;
        i4 = i3+1;
        i5 = (ndx+1)*(ndy+1)+ndx*(i-1)+j; 
        (ndx+1)*i+ndx*(i-1)+j;
        
        a = 4*ndx*(i-1)+4*(j-1)+1;
        b = a+1;
        c = b+1;
        d = c+1;
        element(a,[1 2 3]) = [i5, i1, i2];
        element(b,[1 2 3]) = [i5, i2, i4];
        element(c,[1 2 3]) = [i5, i4, i3];
        element(d,[1 2 3]) = [i5, i3, i1];
    end
end

return