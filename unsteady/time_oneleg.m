function [Vnew,pnew] = time_oneleg(Vn,pn,V_old,p_old,tn,dt,options)

%% one-leg beta method following 
% Symmetry-preserving discretization of turbulent flow, Verstappen and Veldman (JCP 2003)
% or:﻿Direct numerical simulation of turbulence at lower costs (Journal of
% Engineering Mathematics 1997)

% formulation:
% ((beta+1/2)*u^{n+1} -2*beta*u^{n} + (beta-1/2)*u^{n-1})/dt  = 
%  F((1+beta)*u^n - beta*u^{n-1})

%% grid info
% Nu = options.grid.Nu;
% Nv = options.grid.Nv;
Om_inv = options.grid.Om_inv;

%% coefficients of oneleg method
beta = options.time.beta;


%% preprocessing

% store variables at start of time step
% tn     = t;
% Vn     = V;
% pn     = p;
% V = Vn;
% p = pn;

% gradient operator
G      = options.discretization.G;
% divergence operator
M      = options.discretization.M;


%% take time step

% intermediate ('offstep') velocities
t_int  = tn + beta*dt;
V_int  = (1+beta)*Vn - beta*V_old;
p_int  = (1+beta)*pn - beta*p_old; % see paper: 'DNS at lower cost'
%p_temp = p;

% right-hand side of the momentum equation
[~,F_rhs]  = F(V_int,V_int,p_int,t_int,options);

% take a time step with this right-hand side, this gives an 
% intermediate velocity field (not divergence free)
Vtemp = (2*beta*Vn - (beta-0.5)*V_old + dt*Om_inv.*F_rhs)/(beta+0.5);

% to make the velocity field u(n+1) at t(n+1) divergence-free we need
% the boundary conditions at t(n+1)
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn + dt,options);
end
yM = options.discretization.yM;
    

% define an adapted time step; this is only influencing the pressure
% calculation
dt_beta = dt/(beta+0.5);
% divergence of intermediate velocity field is directly calculated with M
f = (M*Vtemp + yM)/dt_beta;

% solve the Poisson equation for the pressure
dp = pressure_poisson(f,tn + dt,options);

% update velocity field
Vnew  = Vtemp - dt_beta*Om_inv.*G*dp;

% update pressure (second order)
pnew  = 2*pn - p_old + (4/3)*dp;

% alternatively, do an additional Poisson solve:
if (options.solversettings.p_add_solve == 1)
    pnew = pressure_additional_solve(Vnew,pn,tn+dt,options);
end




%%
% uhn    = uh;
% vhn    = vh;
% dtn    = dt;
% 
% uh_temp= (1+beta)*uh - beta*uh_old;
% vh_temp= (1+beta)*vh - beta*vh_old;
% 
% cu_temp= uh_temp;
% cv_temp= vh_temp;
% 
% p_temp = (1+beta)*p - beta*p_old; % see paper: 'DNS at lower cost'
% % p_temp = p;
% 
% % convection
% u_ux   = Au_ux*uh_temp+yAu_ux;                 % u at ux
% uf_ux  = Iu_ux*cu_temp+yIu_ux;                 % ubar at ux
% du2dx  = Cux*(uf_ux.*u_ux);    
% 
% u_uy   = Au_uy*uh_temp+yAu_uy;                 % u at uy
% vf_uy  = Iv_uy*cv_temp+yIv_uy;                 % vbar at uy
% duvdy  = Cuy*(vf_uy.*u_uy);
% 
% v_vx   = Av_vx*vh_temp+yAv_vx;                 % v at vx
% uf_vx  = Iu_vx*cu_temp+yIu_vx;                 % ubar at vx
% duvdx  = Cvx*(uf_vx.*v_vx);
% 
% v_vy   = Av_vy*vh_temp+yAv_vy;                 % v at vy
% vf_vy  = Iv_vy*cv_temp+yIv_vy;                 % vbar at vy    
% dv2dy  = Cvy*(vf_vy.*v_vy);
% 
% 
% % diffusion
% d2u    = Diffu*uh_temp + yDiffu;
% d2v    = Diffv*vh_temp + yDiffv;
% 
% 
% % force from old time levels; f( (1+beta)*tn - beta*tn-1, (1+beta)*tn - beta*tn-1) 
% tn=t;
% % 
% t = t+beta*dt;
% force;
% t = tn;
% 
% 
% Ru = 1/(beta+0.5) * ( (2*beta*uh - (beta-0.5)*uh_old) + ...
%                      Omu_inv*dt.*( - du2dx - duvdy + d2u + Fx - Gx*p_temp - y_px ));    
% Rv = 1/(beta+0.5) * ( (2*beta*vh - (beta-0.5)*vh_old) + ...
%                      Omv_inv*dt.*( - duvdx - dv2dy + d2v + Fy - Gy*p_temp - y_py ));      
% 
% R = [Ru;Rv];
% 
% dt = dtn/(beta+0.5);
% 
% pressure_correction;
% 
% dt = dtn;
% 
% % first order pressure:
% % p_old  = p;
% % p      = p + dp;
% 
% % second order pressure:
% p_new  = 2*p - p_old + (4/3)*dp;
% p_old  = p;
% p      = p_new;
% 
% uh_old = uhn;
% vh_old = vhn;
% 
% uh     = V(1:Nu);
% vh     = V(Nu+1:Nu+Nv);