% forward Euler for convection and diffusion
      
% convection
cu = uh;
cv = vh;
convection;

% diffusion
d2u    = Diffu*uh + yDiffu;
d2v    = Diffv*vh + yDiffv;

% force from old time level
force;

Ru = uh - Omu_inv*dt.*( du2dx + duvdy - d2u - Fx + Gx*p + y_px);    
Rv = vh - Omv_inv*dt.*( duvdx + dv2dy - d2v - Fy + Gy*p + y_py);

R = [Ru;Rv];

pressure_correction;

p_old  = p;

uh_old = uh;
vh_old = vh;

uh     = V(1:Nu);
vh     = V(Nu+1:end);
 
p      = p + dp;