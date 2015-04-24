function [EL] = elcook(NX,NY,TYPE)
a = zeros(NY,NX);
for I = 1:NY
    a(I,:) = 1+(I-1)*(NX+1):I*(NX+1)-1;
end
a = sort(reshape(a,size(a,1)*size(a,2),1));
%% ELEMENT MATRIX
% This matrix contains at rows the elements and at columns the nodal points
% relative to each element (in the physical space)
if TYPE == 1 
    EL = zeros(2*NX*NY,3);
    EL(1:2:length(EL),1) = a; 
    EL(1:2:length(EL),2) = a+ones(length(a),1); 
    EL(1:2:length(EL),3) = a+(NX+1)*ones(length(a),1)+ones(length(a),1);
    EL(2:2:length(EL),1) = a;
    EL(2:2:length(EL),2) = a+(NX+1)*ones(length(a),1)+ones(length(a),1);
    EL(2:2:length(EL),3) = EL(2:2:length(EL),2)-ones(length(a),1);
elseif TYPE == 2
    EL = zeros(NX*NY,4);
    EL(1:length(EL),1) = a; 
    EL(1:length(EL),2) = a+ones(length(a),1); 
    EL(1:length(EL),3) = EL(1:length(EL),2)+(NX+1)*ones(length(a),1);
    EL(1:length(EL),4) = a+(NX+1)*ones(length(a),1);
end
end

