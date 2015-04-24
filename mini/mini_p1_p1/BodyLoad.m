function [bf] = BodyLoad(JF,f)

bf(1,1) = f(1,1)*JF/6;
bf(2,1) = f(2,1)*JF/6;
bf(3,1) = f(1,1)*JF/6;
bf(4,1) = f(2,1)*JF/6;
bf(5,1) = f(1,1)*JF/6;
bf(6,1) = f(2,1)*JF/6;
bf(7,1) = f(1,1)*JF;
bf(8,1) = f(2,1)*JF;

return