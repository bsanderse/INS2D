% calculate residual

% convection:
cu = uh;
cv = vh;
convection;

% diffusion
diffusion;

resu   = Omu_inv.*(- convu + d2u - Gx*p - y_px + Fx); 
resv   = Omv_inv.*(- convv + d2v - Gy*p - y_py + Fy);