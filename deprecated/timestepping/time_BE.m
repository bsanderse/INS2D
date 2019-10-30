% implicit convection and diffusion (first order)


% Picard linearization
cu         = uh;
cv         = vh;

% evaluate BC and force at new time 
t = tn + dt;
boundary_conditions;
operator_boundary_conditions;
force;
t = tn;

% convective terms, u-component
% c^n * u^(n+1), c=u
uf_ux      = Iu_ux*cu+yIu_ux;                 % convective velocity, u_bar
C          = Cux*spdiags(uf_ux,0,N1,N1);         
Conv_ux    = C*Au_ux; 
yConv_ux   = C*yAu_ux;

% c^n * u^(n+1), c=v
vf_uy      = Iv_uy*cv+yIv_uy;                 % convective velocity, v_bar
C          = Cuy*spdiags(vf_uy,0,N2,N2);
Conv_uy    = C*Au_uy;
yConv_uy   = C*yAu_uy;

% convective terms, v-component
% c^n * v^(n+1), c=u
uf_vx      = Iu_vx*cu+yIu_vx;                 % convective velocity, u_bar     
C          = Cvx*spdiags(uf_vx,0,N3,N3);
Conv_vx    = C*Av_vx;
yConv_vx   = C*yAv_vx;

% c^n * v^(n+1), c=v
vf_vy      = Iv_vy*cv+yIv_vy;                 % convective velocity, u_bar
C          = Cvy*spdiags(vf_vy,0,N4,N4);       
Conv_vy    = C*Av_vy;
yConv_vy   = C*yAv_vy;      

% force from new time level
tn = t;
t  = t+dt;
force;
t  = tn;

Ru =   ( spdiags(Omu,0,Nu,Nu)/dt + ( - Diffu + Conv_ux + Conv_uy) )...
     \ ( Omu.*uh/dt + ( yDiffu - yConv_ux - yConv_uy + Fx - Gx*p - y_px));       

Rv =   ( spdiags(Omv,0,Nv,Nv)/dt + ( - Diffv + Conv_vx + Conv_vy) ) ...
     \ ( Omv.*vh/dt + ( yDiffv - yConv_vx - yConv_vy + Fy - Gy*p - y_py));
 
R  = [Ru; Rv];

pressure_correction;

p_old  = p;
uh_old = uh;
vh_old = vh;

uh     = V(1:Nu);
vh     = V(Nu+1:end);
p      = p + dp;