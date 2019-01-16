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


nq_u = length(xbi_u);
for qq=1:nq_u

  i = interface_u_q(qq,1);
  j = interface_u_q(qq,2);

  % construct right hand side
  uB = 0.;
  [rhs_u ]   = set_interpolation_rhs(xbi_u,ybi_u,xip_u,yip_u,xin,yp,interface_u_q,u,uB);

  coeff      = squeeze(bilin_mat_inv_u(qq,:,:))*squeeze(rhs_u(qq,:))'; 
  uGC        = 2*uB - coeff'*[1 xip_u(qq) yip_u(qq) xip_u(qq)*yip_u(qq)]';

  f(ind)     = uGC - u(i,j);

end

nq_v = length(xbi_v);
for qq=1:nq_v

  i = interface_v_q(qq,1);
  j = interface_v_q(qq,2);

  % construct right hand side         
  vB = 0.;
  [rhs_v ]   = set_interpolation_rhs(xbi_v,ybi_v,xip_v,yip_v,xp,yin,interface_v_q,v,vB);   

  coeff      = squeeze(bilin_mat_inv_v(qq,:,:))*squeeze(rhs_v(qq,:))';
  vGC        = 2*vB - coeff'*[1 xip_v(qq) yip_v(qq) xip_v(qq)*yip_v(qq)]';

  ind        = Nu + sub2ind(size(v),i,j);
  CD(ind,:)  = 0;
  CD(ind,ind)= 1;

  f(ind)     = vGC - v(i,j);


end
