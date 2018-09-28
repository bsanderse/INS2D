% Adams-Bashforth for convection and Crank-Nicolson for diffusion
      
tn = t;

% evaluate BC and force at starting point for convection
t  = tn;
boundary_conditions;
interpolate_bc;
operator_bc_momentum;
force;

% convection
cu     = uh;
cv     = vh;

convection;

Cu     = du2dx + duvdy;
Cv     = duvdx + dv2dy;

% evaluate BC and force at midpoint for DIFFUSION; note that CN is formally combination
% of F(n) and F(n+1), but this is the same for linear problems 
t = tn + 0.5*dt;
boundary_conditions;
interpolate_bc;
operator_bc_momentum;
force;
t = tn;

% pressure
% p_temp = alfa1*p + alfa2*p_old; % see paper: 'DNS at lower cost'
p_temp = p;


% LU decomposition has been calculated already in
% operator_convection_diffusion

Rur = uh + Omu_inv*dt.*(- (alfa1*Cu + alfa2*Cu_old) + ...
                         + (1-theta)*Diffu*uh + yDiffu + ...
                         + Fx - Gx*p_temp - y_px);

Rvr = vh + Omv_inv*dt.*(- (alfa1*Cv + alfa2*Cv_old) + ...
                         + (1-theta)*Diffv*vh + yDiffv + ...
                         + Fy - Gy*p_temp - y_py);        

if (poisson_diffusion==1)
    b   = L_diffu\Rur;
    Ru  = U_diffu\b;
    
    b   = L_diffv\Rvr;
    Rv  = U_diffv\b;
    
elseif (poisson_diffusion==3)
    [Ru,iter,norm1,norm2]=cg(L_diffu,int64(dia_diffu),int64(ndia_diffu),Rur,...
                            CG_acc,int64(Nu),Ru,int64(CG_maxit));
%                         iter
    [Rv,iter,norm1,norm2]=cg(L_diffv,int64(dia_diffv),int64(ndia_diffv),Rvr,...
                            CG_acc,int64(Nv),Rv,int64(CG_maxit));
%                         iter
end

R = [Ru; Rv];

Cu_old = Cu;
Cv_old = Cv;

% evaluate divergence BC at endpoint, necessary for PC
t = tn + dt;
boundary_conditions;
interpolate_bc;
operator_bc_divergence;
t = tn;

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

pressure_additional_solve;