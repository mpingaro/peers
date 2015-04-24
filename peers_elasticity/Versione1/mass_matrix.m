%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     MASS MATRIX OF PEERS ELEMENT                        %
%  by Dott. Ing. Marco Pingaro (PhD student)                              %
%  Department of Civil Engineering and Architecture,                      %
%             Via Ferrata 3, 27100, Pavia                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  Computing the mass matrix of element PEERS                             %
%                                                                         %
%  INPUT : P is a vector 6x1                                              %
%          rho is a density per unit area of the material                 %
%  OUTPUT: mass is a elementar mass matrix                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mass = mass_matrix(P,rho)

% Area 
A = 0.5*(( P(3,1)-P(1,1) )*( P(6,1)-P(2,1) )...
    -( P(5,1)-P(1,1) )*( P(4,1)-P(2,1) ));

mass = sparse(2,2);
mass(1,1) = rho*A;
mass(2,2) = mass(1,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%