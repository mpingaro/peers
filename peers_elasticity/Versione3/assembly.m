function [K,p] = assembly(coordinates,element,mc,mu,lambda,f)

ngdls = max(max(mc));
nelem = size(element,1);
nnod = size(coordinates,1);

% Preallocazione matrici e vettori
bf = sparse(2*nelem,1);
bg = sparse(ngdls,1);
bh = sparse(nnod,1);
A = sparse(ngdls,ngdls); 
B = sparse(ngdls,2*nelem);
C = sparse(ngdls,nnod);

[psi,eta] = EvaluateShape;
% ASSEMBLY GLOBAL MATRIX A B & C
for k=1:nelem
    P([1 2])=coordinates(element(k,1),[1 2]);
    P([3 4])=coordinates(element(k,2),[1 2]);
    P([5 6])=coordinates(element(k,3),[1 2]);
    s = mc(k,9);
   
    [AREA,DF,DFF,JF] = TriTrasformation(P);
    [psiF,etaF] = EvaluateShapePhysics(DF,DFF,JF,psi,eta);
    
    AELEM = SigmaTauNEW(psiF,JF,mu,lambda,s);
    BELEM = UdivTau(s);
    CELEM = EtaTauNEW(psiF,etaF,JF,s);
    
    b = bodyload(AREA,f);
    bf(2*(k-1)+[1 2]) = b([1 2],1);
    for i=1:8
        for j=1:8
            if i>j
                AELEM(i,j)=AELEM(j,i);
            end
            A(mc(k,i),mc(k,j)) = A(mc(k,i),mc(k,j))+AELEM(i,j);
        end
        B(mc(k,i),2*(k-1)+[1 2]) = BELEM(i,[1 2]);
        for z =1:3
            C(mc(k,i), element(k,z)) = C(mc(k,i),element(k,z))...
                +CELEM(i,z);
        end
    end
end
% ASSEMBLY GLOBAL SISTEM K & GLOBAL LOAD VECTOR P
K = [A, B, C; B', sparse(2*nelem,2*nelem+nnod); C', sparse(nnod,2*nelem+nnod)];
p = [bg;bf;bh];
end
