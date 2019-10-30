% implicit convection and diffusion (second order)
% almost same as time_IM1, except for cu and cv at n and n+1

CN         = 0.5;

cu         = uh;
cv         = vh;

% pressure
% p_temp     = p;
p_temp     = 1.5*p - 0.5*p_old;

%% convective and diffusive terms from old time level
Conv_u     = Cux*((Iu_ux*cu+yIu_ux).*(Au_ux*uh+yAu_ux)) + ...
             Cuy*((Iv_uy*cv+yIv_uy).*(Au_uy*uh+yAu_uy));
         
Conv_v     = Cvx*((Iu_vx*cu+yIu_vx).*(Av_vx*vh+yAv_vx)) + ...
             Cvy*((Iv_vy*cv+yIv_vy).*(Av_vy*vh+yAv_vy));
         
Diff_u     = Diffu*uh;
Diff_v     = Diffv*vh;

%% convective and diffusive terms at new time level
% Picard (first order):
% cu         = uh;
% cv         = vh;

% extrapolated Picard (second order) for n+1:
cu         = 2*uh - uh_old;
cv         = 2*vh - vh_old;

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


% force is a combination from old and new time level
tn = t;
force; % time level n
Fx_old = Fx; 
Fy_old = Fy; 
t  = t+dt;
force; % time level n+1 (if the force is not known we need to linearize!)
Fx = 0.5*(Fx+Fx_old);
Fy = 0.5*(Fy+Fy_old);
t  = tn;


Ru =   (  spdiags(Omu,0,Nu,Nu)/dt + CN*( - Diffu + Conv_ux + Conv_uy) )...
     \ (  Omu.*uh/dt + CN*(Diff_u - Conv_u - yConv_ux - yConv_uy) + ...
                          (yDiffu + Fx - Gx*p_temp - y_px));       

Rv =   (  spdiags(Omv,0,Nv,Nv)/dt + CN*( - Diffv + Conv_vx + Conv_vy) ) ...
     \ (  Omv.*vh/dt + CN*(Diff_v - Conv_v - yConv_vx - yConv_vy) + ...
                          (yDiffv + Fy - Gy*p_temp - y_py));
 
R = [Ru; Rv];

pressure_correction;

% first order pressure:
% p_old  = p;
% p      = p + dp;

% second order pressure:
p_new  = 2*p - p_old + (4/3)*dp;
p_old  = p;
p      = p_new;

uh_old = uh;
vh_old = vh;

uh     = V(1:Nu);
vh     = V(Nu+1:Nu+Nv);