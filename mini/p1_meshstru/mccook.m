function [MCORR] = mccook(EL,TYPE)
%% CORRESPONDENCE MATRIX
% In the correspondence matrix rows are relative to elements whereas 
% columns are relative to the degrees of freedom for each element.
if TYPE == 1
    i = EL(:,1)'; j = EL(:,2)'; z = EL(:,3)';
    MCORR = [2*i-1 ; 2*i ; 2*j-1 ; 2*j ; 2*z-1 ; 2*z]; MCORR = MCORR';
elseif TYPE == 2
    i = EL(:,1)'; j = EL(:,2)'; t = EL(:,3)'; z = EL(:,4)';
    MCORR = [2*i-1 ; 2*i ; 2*j-1 ; 2*j ; 2*t-1 ; 2*t ; 2*z-1 ; 2*z]; MCORR = MCORR';
end
end

