%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            VALUE OF STRESS                              %       
%  by Dott. Ing. Marco Pingaro (PhD student)                              %
%  Department of Civil Engineering and Architecture,                      %
%             Via Ferrata 3, 27100, Pavia                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  ValStr = StressValue(mc,st,nel)

ValStr = zeros(3*nel,4);
for i =1:nel
    s([1 2 3 4 5 6 7 8],1) = st(mc(i,[1 2 3 4 5 6 7 8]),1);
    x = [0 1 0; 0 0 1];
    for j = 1:3
        x(1) = x(1,j);
        x(2) = x(2,j);
        psi(:,1) = [x(1); -1+x(2)];
        psi(:,2) = [(2/sqrt(2))*x(1); (2/sqrt(2))*x(2)];
        psi(:,3) = [-1+x(1); x(2)];
        % BOLLA
        dBx = (x(2)-x(2)*x(2)-2*x(1)*x(2))*27; 
        dBy = (x(1)-x(1)*x(1)-2*x(1)*x(2))*27;
        psi(:,4) = [dBy; -dBx];
        % SIGMA XX
        ValStr(3*(i-1)+j,1) = s(1,1)*psi(1,1)+ s(3,1)*psi(1,2)+...
            s(5,1)*psi(1,3)+s(7,1)*psi(1,4); 
        % SIGMA XY
        ValStr(3*(i-1)+j,2) = s(1,1)*psi(2,1)+ s(3,1)*psi(2,2)+...
            s(5,1)*psi(2,3)+s(7,1)*psi(2,4);
        % SIGMA YX
        ValStr(3*(i-1)+j,3) = s(2,1)*psi(1,1)+ s(4,1)*psi(1,2)+...
            s(6,1)*psi(1,3)+s(8,1)*psi(1,4);
        % SIGMA YY
        ValStr(3*(i-1)+j,4) = s(2,1)*psi(2,1)+ s(4,1)*psi(2,2)+...
            s(6,1)*psi(2,3)+s(8,1)*psi(2,4);
    end    
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%