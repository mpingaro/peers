function plt_solution(msh,solu)
%
% Contiene tutti i plottaggi e i salvataggi
%%
figure
plotmesh(msh.element,msh.coordinates)
title('Undeformed Mesh')
%print('undeformed_msh.eps')
%
figure
plotmesh(msh.element,solu.defo)
title('Deformed Mesh')
%print('deformed_msh.eps') 
%
figure
ShowDisplacement(msh.element,solu.defo,solu.ux)
title('Displacement U_{x}')
%print('displacement_ux.eps')
%
figure
ShowDisplacement(msh.element,solu.defo,solu.uy)
title('Displacement U_{y}')
%print('displacement_uy.eps')
%
figure
ShowRotRigid(msh.element,msh.coordinates,solu.etasp)
title('Displacement Rotation')
%print('rotation.eps')
%
figure
ShowStress(msh.element,msh.coordinates,solu.ValStr(:,1));
title('Stress Sigma_{xx}')
%print('stress_xx.eps')
%
figure
ShowStress(msh.element,msh.coordinates,solu.ValStr(:,4));
title('Stress Sigma_{yy}')
%print('stress_yy.eps')
%
figure
ShowStress(msh.element,msh.coordinates,solu.ValStr(:,2));
title('Stress Sigma_{xy}')
%print('stress_xy.eps')
%
figure
ShowStress(msh.element,msh.coordinates,solu.ValStr(:,3));
title('Stress Sigma_{yx}')
%print('stress_yx.eps')

end
%% END FUNCTION