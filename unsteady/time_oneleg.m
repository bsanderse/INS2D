% one-leg beta method
uhn    = uh;
vhn    = vh;
dtn    = dt;

uh_temp= (1+beta)*uh - beta*uh_old;
vh_temp= (1+beta)*vh - beta*vh_old;

cu_temp= uh_temp;
cv_temp= vh_temp;

p_temp = (1+beta)*p - beta*p_old; % see paper: 'DNS at lower cost'
% p_temp = p;

% convection
u_ux   = Au_ux*uh_temp+yAu_ux;                 % u at ux
uf_ux  = Iu_ux*cu_temp+yIu_ux;                 % ubar at ux
du2dx  = Cux*(uf_ux.*u_ux);    

u_uy   = Au_uy*uh_temp+yAu_uy;                 % u at uy
vf_uy  = Iv_uy*cv_temp+yIv_uy;                 % vbar at uy
duvdy  = Cuy*(vf_uy.*u_uy);

v_vx   = Av_vx*vh_temp+yAv_vx;                 % v at vx
uf_vx  = Iu_vx*cu_temp+yIu_vx;                 % ubar at vx
duvdx  = Cvx*(uf_vx.*v_vx);

v_vy   = Av_vy*vh_temp+yAv_vy;                 % v at vy
vf_vy  = Iv_vy*cv_temp+yIv_vy;                 % vbar at vy    
dv2dy  = Cvy*(vf_vy.*v_vy);


% diffusion
d2u    = Diffu*uh_temp + yDiffu;
d2v    = Diffv*vh_temp + yDiffv;


% force from old time levels; f( (1+beta)*tn - beta*tn-1, (1+beta)*tn - beta*tn-1) 
tn=t;
% 
t = t+beta*dt;
force;
t = tn;


Ru = 1/(beta+0.5) * ( (2*beta*uh - (beta-0.5)*uh_old) + ...
                     Omu_inv*dt.*( - du2dx - duvdy + d2u + Fx - Gx*p_temp - y_px ));    
Rv = 1/(beta+0.5) * ( (2*beta*vh - (beta-0.5)*vh_old) + ...
                     Omv_inv*dt.*( - duvdx - dv2dy + d2v + Fy - Gy*p_temp - y_py ));      

R = [Ru;Rv];

dt = dtn/(beta+0.5);

pressure_correction;

dt = dtn;

% first order pressure:
% p_old  = p;
% p      = p + dp;

% second order pressure:
p_new  = 2*p - p_old + (4/3)*dp;
p_old  = p;
p      = p_new;

uh_old = uhn;
vh_old = vhn;

uh     = V(1:Nu);
vh     = V(Nu+1:Nu+Nv);