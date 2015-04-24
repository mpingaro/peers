function msh = GeoTri(point)
% -- GENERA STRUTTURA DATI GEOMETRIA TRIANGOLO ---------------------------%
% -- by Dott. Ing. Marco Pingaro (PhD Students IUSS) 
% -- Created in 2013
%
% -- CALCOLA DATI I PUNTI DEL TRIANGOLO:
%       BARICENTRO
%       AREA
%       MATRICE JACOBIANA
%       DETERMINANTE DELLO JACOBIANO

% --- Baricentro Triangolo
msh.bari(1) = sum(point(:,1))/3;
msh.bari(2) = sum(point(:,2))/3;

% --- Matrice Jacobiana
msh.jac(1,1) = point(2,1)-point(1,1);
msh.jac(1,2) = point(3,1)-point(1,1);
msh.jac(2,1) = point(2,2)-point(1,2);
msh.jac(2,2) = point(3,2)-point(1,2);

% --- Determinante della matrice Jacobiana
msh.djac = msh.jac(1,1)*msh.jac(2,2)-msh.jac(1,2)*msh.jac(2,1);

% --- Area del triangolo
msh.area = msh.djac/2;

end
% -- END FUNCTION