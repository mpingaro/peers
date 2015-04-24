function [gauss, weight] = Gauss_Quadrature(order)

if order == 1
    %% QUADRATURA DI GAUSS (grado di precisione 1, 1 points)
    % Weight of quadratura
    weight(1,1) = 0.500000000000000;
    % Point value
    gauss(1,1) = 0.5;
    gauss(2,1) = 0.5;
elseif order == 2
    %% QUADRATURA DI GAUSS (grado di precisione 2, 3 points)
    % Weight of quadrature
    weight(1,1) = 0.166666666666667;
    weight(1,2) = 0.166666666666667;
    weight(1,3) = 0.166666666666667;
    % Point value
    gauss(1,1) = 0.5;
    gauss(2,1) = 0.0;
    gauss(1,2) = 0.5;
    gauss(2,2) = 0.5;
    gauss(1,3) = 0.0;
    gauss(2,3) = 0.5;
elseif order == 5 
    %% QUADRATURA DI GAUSS (grado di precisione 5, 7 points)
    % Weight of quadrature
    weight(1,1) = 0.112500000000000; weight(1,2) = 0.066197076394253; 
    weight(1,3) = 0.066197076394253; weight(1,4) = 0.066197076394253; 
    weight(1,5) = 0.062969590272414; weight(1,6) = 0.062969590272414; 
    weight(1,7) = 0.062969590272414;
    % Point value 
    gauss(1,1) = 0.333333333333333; gauss(1,2) = 0.470142064105115; 
    gauss(1,3) = 0.059715871789770; gauss(1,4) = 0.470142064105115; 
    gauss(1,5) = 0.101286507323456; gauss(1,6) = 0.797426985353087;
    gauss(1,7) = 0.101286507323456; 
    gauss(2,1) = 0.333333333333333; gauss(2,2) = 0.470142064105115; 
    gauss(2,3) = 0.470142064105115; gauss(2,4) = 0.059715871789770;
    gauss(2,5) = 0.101286507323456; gauss(2,6) = 0.101286507323456;
    gauss(2,7) = 0.797426985353087;
end


end