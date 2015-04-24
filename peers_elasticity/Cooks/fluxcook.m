function [FLUX] = fluxcook(NX,NY)

FLUX = zeros(2*NX*NY,8);
I8 = 6*NX*NY+2*(NX+NY)+1;  % index g.d.l bubble
for I = 1:NY
    for J = 1:NX
        I1 = (I-1)*(2*NX+NX+1)+J; 
        I4 = I1+3*NX+1; I5 = I1+NX-1+J;
        I6 = I1+NX+J;
        FLUX(1+(I-1)*2*NX+(J-1)*2,1) = I1;
        FLUX(1+(I-1)*2*NX+(J-1)*2,2) = I1+NX+1+J;
        FLUX(1+(I-1)*2*NX+(J-1)*2,3) = I1+NX+J;
        FLUX(2+(I-1)*2*NX+(J-1)*2,[1 2 3 4 5 6 7 8]) =...
            [I4*2-1, I4*2, I5*2, I5*2, I6*2-1, I6*2, I8+2 , I8+3];
        I8 = I8 + 4;
    end
end

end

