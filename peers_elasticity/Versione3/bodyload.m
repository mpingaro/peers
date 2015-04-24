function     [b] = bodyload(A,f)

b(1,1) = f(1,1)*A;
b(2,1) = f(2,1)*A;

end