function [p] = Load(f,ngdl,nelem);

% ASSEMBLY VECTOR p (LOAD VECTOR)
bg = sparse(ngdl,1);       % Dirichlet condition
bf = sparse(2*nelem,1);    % Body forces
bh = sparse(3*nelem,1);    % Simmetry imposition 

for k=1:nelem
    bf(2*(k-1)+[1 2],1) =[f(1,1), f(2,1)];
end

p = [bg;bf;bh];

return