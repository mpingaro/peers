function [psi,eta] = EvaluateShape

% Gauss Evaluetion Points
q1=1/3; q2=(6+sqrt(15))/21; q3=(9-2*sqrt(15))/21;
q4=(6-sqrt(15))/21; q5=(9+2*sqrt(15))/21;
xnod=[q1 q2 q3 q2 q4 q5 q4]; 
ynod=[q1 q2 q2 q3 q4 q4 q5];

for i=1:7
    % Stress
    psi(1,i) = xnod(i); 
    psi(2,i) = -1+ynod(i);
    psi(3,i) = 2/sqrt(2)*xnod(i);
    psi(4,i) = 2/sqrt(2)*ynod(i);
    psi(5,i) = -1+xnod(i);
    psi(6,i) = ynod(i);
    psi(7,i) = (ynod(i)-ynod(i)*ynod(i)-2*xnod(i)*ynod(i))*120;
    psi(8,i) = (xnod(i)-xnod(i)*xnod(i)-2*xnod(i)*ynod(i))*120;
    % Rotations
    eta(1,i) = 1-xnod(i)-ynod(i);
    eta(2,i) = xnod(i);
    eta(3,i) = ynod(i);
    
end

return