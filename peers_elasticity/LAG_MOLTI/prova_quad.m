% INTEGRALE SUL TRIANGOLO PARENTE

I =0;
for k=1:7
    [gauss, weight] = Gauss_Quadrature(5);
    x(1) = gauss(1,k); x(2) = gauss(2,k);
    
    I = I + weight(k)*x(1)*x(2)*( 1-x(1)-x(2) );
end


disp(I)