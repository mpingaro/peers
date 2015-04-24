%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    COORDINATES OF DEFORMED MESH                         %       
%  by Dott. Ing. Marco Pingaro (PhD student)                              %
%  Department of Civil Engineering and Architecture,                      %
%             Via Ferrata 3, 27100, Pavia                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function defo = DefoMesh(coord, ux, uy, nodes2element)

nnd = size(coord,1);
defo = zeros(nnd,2);

for i = 1:nnd
    spo_x = 0;
    spo_y = 0;
    ndiv = 1;
    for ii = 1:nnd
        if nodes2element(i,ii) ~= 0
            nel = nodes2element(i,ii);
            ndiv = ndiv + 1;
            spo_x = spo_x + ux(nel,1);
            spo_y = spo_y + uy(nel,1);
        end
    end
    ndiv = ndiv-1;
    spo_x = spo_x/ndiv;
    spo_y = spo_y/ndiv;
    defo(i,1) = coord(i,1) + spo_x;
    defo(i,2) = coord(i,2) + spo_y;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%