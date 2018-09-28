%% regularization matrices

alfa     = (1/16)*deltax^2;

% diffusive matrices in finite-difference setting, without viscosity
Diffu_f  = spdiags(Omu_inv,0,Nu,Nu)*(Dux*Su_ux + Duy*Su_uy);
Diffv_f  = spdiags(Omv_inv,0,Nv,Nv)*(Dux*Sv_vx + Dvy*Sv_vy);

yDiffu_f = Omu_inv.*(Dux*ySu_ux + Duy* ySu_uy);
yDiffv_f = Omv_inv.*(Dvx*ySv_vx + Dvy* ySv_vy);
