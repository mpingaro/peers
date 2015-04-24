function u = Solve_PEERS(geo,alpha,beta,bf)

%% ASSEMBLY GLOBAL MATRIX
B = sparse(8*geo.nelem,8*geo.nelem);
C = sparse(8*geo.nelem,2*geo.nelem);
D = sparse(8*geo.nelem,geo.nnod);
E = sparse(8*geo.nelem,2*size(geo.intedge,1));

ngdlt = 10*geo.nelem+geo.nnod+2*size(geo.intedge,1);
load = sparse(ngdlt,1);

for k = 1:geo.nelem
    point(1,[1 2]) = geo.coordinates( geo.element(k,1),[1 2] );
    point(2,[1 2]) = geo.coordinates( geo.element(k,2),[1 2] );
    point(3,[1 2]) = geo.coordinates( geo.element(k,3),[1 2] );
    msh = GeoTri(point);
    [BELEM,CELEM,DELEM] = PEERS_Element(msh,point,alpha,beta);
    for i=1:8
        for j=1:8
            B(8*(k-1)+i,8*(k-1)+j) = BELEM(i,j);
        end
        C(8*(k-1)+i,2*(k-1)+[1 2]) = CELEM(i,[1 2]);
        for ii = 1:3
            D(8*(k-1)+i, geo.element(k,ii) ) = D(8*(k-1)+i, geo.element(k,ii)) + DELEM(i,ii);
        end
    end
    load(8*geo.nelem+2*(k-1)+[1 2],1) = Bulk_load(msh,bf);
end

% Continuità lati interni
for k = 1: size(geo.intedge,1)
    [EELEM] = PEERS_intcontinuity(geo,geo.intedge(k));
    for i = 1:8
        E(8*(geo.edg2elm(geo.intedge(k),3)-1)+i, 2*(k-1)+[1 2]) = EELEM(i,[1 2],1);
        E(8*(geo.edg2elm(geo.intedge(k),4)-1)+i, 2*(k-1)+[1 2]) = EELEM(i,[1 2],2);
    end
end

% Assemblo sistema globale
stiff = [B, C, D, E; 
    C', sparse(2*geo.nelem, 2*geo.nelem+geo.nnod+2*size(geo.intedge,1));
    D', sparse(geo.nnod, 2*geo.nelem+geo.nnod+2*size(geo.intedge,1));
    E', sparse(2*size(geo.intedge,1), 2*geo.nelem+geo.nnod+2*size(geo.intedge,1))];

% Solve
soluz = stiff\load;

% Estrazione risultati
u.sigma = soluz(1:8*geo.nelem); 
u.disp = reshape( soluz(8*geo.nelem+1:10*geo.nelem), 2, []);
u.rot = soluz(10*geo.nelem+1:10*geo.nelem+geo.nnod);
u.molt = soluz(10*geo.nelem+geo.nnod+1:10*geo.nelem+geo.nnod+2*size(geo.intedge,1));

end