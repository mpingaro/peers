function [FLUX] = fluxcook(NX,NY,TYPE)

if TYPE == 1
    FLUX = zeros(2*NX*NY,4);
    for I = 1:NY
        for J = 1:NX
            I1 = (I-1)*(2*NX+NX+1)+J; 
            I4 = I1+3*NX+1; I5 = I1+NX-1+J;
            I6 = I1+NX+J;
            FLUX(1+(I-1)*2*NX+(J-1)*2,1) = I1;
            FLUX(1+(I-1)*2*NX+(J-1)*2,2) = I1+NX+1+J;
            FLUX(1+(I-1)*2*NX+(J-1)*2,3) = I1+NX+J;
            FLUX(2+(I-1)*2*NX+(J-1)*2,[1 2 3]) = [I4 I5 I6];            
        end
    end
elseif TYPE == 2
    FLUX = zeros(NX*NY,5);
    for I = 1:NY
        for J = 1:NX
            I1 = (I-1)*(2*NX+1)+J;
            I2 = (I-1)*(2*NX+1)+NX+1+J;
            I3 = I*(2*NX+1)+J; I4 = I2-1;
            FLUX((I-1)*NX+J,[1 2 3 4]) = [I1 I2 I3 I4];
        end
    end
end
end

