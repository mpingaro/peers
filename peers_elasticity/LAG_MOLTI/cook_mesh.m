function geo = cook_mesh(sizemesh)
% Generazione mesh non trutturata MENSOLA DI COOKs.
pv =[0.0, 0.0; 48.0, 44.0; 48.0, 60.0; 0.0, 44.0; 0.0, 0.0];
[coordinates,element] = ...
    distmesh2d(@dpoly,@huniform,sizemesh,[0,0; 48,60],pv,pv);

geo.coordinates = coordinates;
geo.element = element;
geo.nelem = size(element,1);
geo.nnod = size(coordinates,1);

end