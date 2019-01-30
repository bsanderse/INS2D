%% Adams-Bashforth for convection and Crank-Nicolson for diffusion
% formulation:
% (u^{n+1} - u^{n})/dt  = -(alfa1*(conv^n) + alfa2*(conv^{n-1})) +
%                           theta*diff^{n+1} + (1-theta)*diff^{n} +
%                           theta*F^{n+1}    + (1-theta)*F^{n}
%                           theta*bc^{n+1}   + (1-theta)*bc^{n} 
%                           - G*p + y_p
% where bc are boundary conditions of diffusion

% rewrite as:
% (I/dt - theta*D)*u^{n+1} = (I/dt - (1-theta)*D)*u^{n} + 
%                           -(alfa1*(conv^n) + alfa2*(conv^{n-1})) +
%                            theta*F^{n+1}    + (1-theta)*F^{n}
%                            theta*bc^{n+1} + (1-theta)*bc^{n}
%                           - G*p + y_p

% the LU decomposition of the first matrix is precomputed in
% operator_convection_diffusion

%%
% Adams-Bashforth coefficients
alfa1 = 3/2;
alfa2 = -1/2;

% CN coefficients
% theta = 1/2;
theta  = options.time.theta;
% store variables at start of time step
tn     = t;
Vn     = V;
pn     = p;
uh     = V(1:Nu);
vh     = V(Nu+1:end);

% grid variables
Omu_inv = options.grid.Omu_inv;
Omv_inv = options.grid.Omv_inv;

%% evaluate BC and force at starting point
[Fx1,Fy1]      = force(tn,options);
% unsteady BC at current time
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn,options);
end

% convection of current solution
[convu, convv] = convection(V,tn,options,0);

% diffusion of current solution
Diffu  = options.discretization.Diffu;
Diffv  = options.discretization.Diffv;
yDiffu1 = options.discretization.yDiffu;
yDiffv1 = options.discretization.yDiffv;

%% evaluate BC and force at end of time step

% unsteady BC at next time
[Fx2,Fy2]      = force(tn+dt,options);
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn+dt,options);
end
% diffusion BC at new time level
yDiffu2 = options.discretization.yDiffu;
yDiffv2 = options.discretization.yDiffv;

Gx   = options.discretization.Gx;
Gy   = options.discretization.Gy;
y_px = options.discretization.y_px;
y_py = options.discretization.y_py;


%% Crank-Nicolson weighting for force and diffusion boundary conditions
Fx = (1-theta)*Fx1 + theta*Fx2;
Fy = (1-theta)*Fy1 + theta*Fy2;
yDiffu = (1-theta)*yDiffu1 + theta*yDiffu2;
yDiffv = (1-theta)*yDiffv1 + theta*yDiffv2;


% pressure
% p_temp = alfa1*p + alfa2*p_old; % see paper: 'DNS at lower cost'
% p_temp = p;

%% right hand side of the momentum equation update
Rur = uh + Omu_inv*dt.*(- (alfa1*convu + alfa2*convu_old) + ...
                         + (1-theta)*Diffu*uh + yDiffu + ...
                         + Fx - Gx*p - y_px);

Rvr = vh + Omv_inv*dt.*(- (alfa1*convv + alfa2*convv_old) + ...
                         + (1-theta)*Diffv*vh + yDiffv + ...
                         + Fy - Gy*p - y_py);        

% LU decomposition of diffusion part has been calculated already in
% operator_convection_diffusion

if (poisson_diffusion==1)
    b   = options.discretization.L_diffu\Rur;
    Ru  = options.discretization.U_diffu\b;
    
    b   = options.discretization.L_diffv\Rvr;
    Rv  = options.discretization.U_diffv\b;
    
elseif (poisson_diffusion==3)
    [Ru,iter,norm1,norm2]=cg(L_diffu,int64(dia_diffu),int64(ndia_diffu),Rur,...
                            CG_acc,int64(Nu),Ru,int64(CG_maxit));
%                         iter
    [Rv,iter,norm1,norm2]=cg(L_diffv,int64(dia_diffv),int64(ndia_diffv),Rvr,...
                            CG_acc,int64(Nv),Rv,int64(CG_maxit));
%                         iter
end

R = [Ru; Rv];

convu_old = convu;
convv_old = convv;

% evaluate divergence BC at endpoint, necessary for PC
yM      = options.discretization.yM;

pressure_correction;

% first order pressure:
% p_old  = p;
p      = p + dp;

% second order pressure:
% p_new  = 2*p - p_old + (4/3)*dp;
% p_old  = p;
% p      = p_new;

if (options.solversettings.p_add_solve == 1)
    p = pressure_additional_solve(V,p,tn+dt,options);
end