function B = UdivTau(s)
% -- Matrix u*div(Tau)
B=sparse(8,2);
%
B(1,1) = s;
B(2,2) = s;
B(3,1) = 2/sqrt(2)*s;
B(4,2) = 2/sqrt(2)*s;
B(5,1) = s;
B(6,2) = s;

return