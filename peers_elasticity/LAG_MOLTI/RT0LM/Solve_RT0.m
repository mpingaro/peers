function u = Solve_RT0(geo,bf)

%% ASSEMBLY GLOBAL MATRIX
B = sparse(3*geo.nelem,3*geo.nelem);
C = sparse(3*geo.nelem,geo.nelem);
D = sparse(3*geo.nelem,size(geo.intedge,1));

ngdlt = 4*geo.nelem+size(geo.intedge,1);
load = sparse(ngdlt,1);

for k = 1:geo.nelem
    point(1,[1 2]) = geo.coordinates( geo.element(k,1),[1 2] );
    point(2,[1 2]) = geo.coordinates( geo.element(k,2),[1 2] );
    point(3,[1 2]) = geo.coordinates( geo.element(k,3),[1 2] );
    msh = GeoTri(point);
    [BELEM,CELEM] = RT0_Element(msh,point);
    for i=1:3
        for j=1:3
            B(3*(k-1)+i,3*(k-1)+j) = BELEM(i,j);
        end
        C(3*(k-1)+i, k) = CELEM(i,1);
    end
    load(3*geo.nelem + k,1) = -msh.area*bf;
end

for k = 1: size(geo.intedge,1)
    [DELEM] = RT0_intcontinuity(geo,geo.intedge(k));
    for i = 1:3
        D(3*(geo.edg2elm(geo.intedge(k),3)-1)+i, k) = DELEM(i,1,1);
        D(3*(geo.edg2elm(geo.intedge(k),4)-1)+i, k) = DELEM(i,1,2);
    end
end
% Assemblo sistema globale
stiff = [B, C, D; 
    C', sparse(geo.nelem, geo.nelem+size(geo.intedge,1));
    D', sparse(size(geo.intedge,1), geo.nelem+size(geo.intedge,1))];

% Solve
soluz = stiff\load;

% Estrazione risultati
u.sigma = soluz(1:3*geo.nelem); 
u.disp = soluz(3*geo.nelem+1:4*geo.nelem);
u.molt = soluz(4*geo.nelem+1:4*geo.nelem+size(geo.intedge,1));

end