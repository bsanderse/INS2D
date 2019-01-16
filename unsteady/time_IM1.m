% Implicit Midpoint without iteration to remove linearization error
% pressure correction

CN     = 0.5;

% cu and cv should be the same for time n and time n+1 to conserve energy
% time_IM1 and time_CN are almost equal, but with cu, cv different for n and n+1

% extrapolated Picard:
cu         = 1.5*uh - 0.5*uh_old;
cv         = 1.5*vh - 0.5*vh_old; 

% pressure
p_temp     = 1.5*p - 0.5*p_old;
    

% evaluate BC and force at intermediate time 
t = tn + 0.5*dt;
boundary_conditions;
operator_boundary_conditions;
force;
t = tn;
    
%% convective and diffusive terms from old time level
Conv_u     = Cux*((Iu_ux*cu+yIu_ux).*(Au_ux*uh+yAu_ux)) + ...
             Cuy*((Iv_uy*cv+yIv_uy).*(Au_uy*uh+yAu_uy));

Conv_v     = Cvx*((Iu_vx*cu+yIu_vx).*(Av_vx*vh+yAv_vx)) + ...
             Cvy*((Iv_vy*cv+yIv_vy).*(Av_vy*vh+yAv_vy));

Diff_u     = Diffu*uh;
Diff_v     = Diffv*vh;
         

%% convective terms, u-component
% c^n * u^(n+1), c=u
uf_ux      = Iu_ux*cu+yIu_ux;                     % convective velocity, u_bar
C          = Cux*spdiags(uf_ux,0,N1,N1);         
Conv_ux    = C*Au_ux; 
yConv_ux   = C*yAu_ux;

% c^n * u^(n+1), c=v
vf_uy      = Iv_uy*cv+yIv_uy;                     % convective velocity, v_bar
C          = Cuy*spdiags(vf_uy,0,N2,N2);
Conv_uy    = C*Au_uy;
yConv_uy   = C*yAu_uy;


%% convective terms, v-component
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

% force is in between time levels
tn = t;
t  = t+0.5*dt;
force; % time level n
t  = tn;


Ru =  (  spdiags(Omu,0,Nu,Nu)/dt + CN*( - Diffu + Conv_ux + Conv_uy) )...
    \ (  Omu.*uh/dt + CN*(Diff_u - Conv_u - yConv_ux - yConv_uy) + ...
                         (yDiffu  + Fx - Gx*p_temp - y_px));       

Rv =  (  spdiags(Omv,0,Nv,Nv)/dt + CN*( - Diffv + Conv_vx + Conv_vy) ) ...
    \ (  Omv.*vh/dt + CN*(Diff_v - Conv_v - yConv_vx - yConv_vy ) + ...
                         (yDiffv + Fy - Gy*p - y_py));

R = [Ru; Rv];
    
pressure_correction;

% first order pressure:
p_old  = p;
p      = p + dp;

% second order pressure:
% p_new  = 2*p - p_old + (4/3)*dp;
% p_old  = p;
% p      = p_new;

uh_old = uh;
vh_old = vh;

uh     = V(1:Nu);
vh     = V(Nu+1:Nu+Nv);