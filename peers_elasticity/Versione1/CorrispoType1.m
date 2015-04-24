function [mc]=CorrispoType1(ndx,ndy)

nelem = 2*ndx*ndy;           % numero elementi
mc = zeros(nelem,9);         % prealloco matrice di corrispondenza
i6 = 3*nelem+2*(ndx+ndy)+1;  % indice g.d.l bolla
for i=1:ndy
    for j = 1:ndx
        i1 = (6*ndx+2)*(i-1)+2*(j-1)+1;
        i2 = (6*ndx+2)*(i-1)+4*(j-1)+ndx*2+1;
        i3 = i2+2;
        i4 = i2+4;
        i5 = (6*ndx+2)*(i-1)+2*(j-1)+6*ndx+3;
        
        mc(2*ndx*(i-1)+2*(j-1)+1,[1 2 3 4 5 6 7 8]) =...
            [i1, i1+1, i3, i3+1, i2, i2+1, i6, i6+1];
        mc(2*ndx*(i-1)+2*(j-1)+2,[1 2 3 4 5 6 7 8]) =...
            [i5, i5+1, i3, i3+1, i4, i4+1, i6+2, i6+3];
        i6 = i6+4;
    end
end
% indice dei segni flussi
segno = zeros(nelem,1);
segno([1:2:nelem],1) = 1; segno([2:2:nelem],1) = -1;
mc(:,9) = segno(:);

%save mc
return