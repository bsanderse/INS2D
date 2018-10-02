% general explicit Runge-Kutta method
% integration of boundary values

Nu = options.grid.Nu;
Nv = options.grid.Nv;
Np = options.grid.Np;

Omu_inv = options.grid.Omu_inv;
Omv_inv = options.grid.Omv_inv;

Gx = options.discretization.Gx;
Gy = options.discretization.Gy;


% get coefficients of RK method
if (isnumeric(options.time.RK))
    options.time.RK = num2str(options.time.RK);
end
[A_RK,b_RK,c_RK] = getRKmethod(options.time.RK);
% RK_order = check_orderconditions(A_RK,b_RK,c_RK);

% number of stages
s_RK = length(b_RK);

% we work with the following 'shifted' Butcher tableau, because A_RK(1,1)
% is always zero for explicit methods
A_RK = [A_RK(2:end,:); b_RK'];

% number of stages
c_RK = [c_RK(2:end);1]; % 1 is the time level of final step

% store variables at start of time step
tn     = t;
uhn    = uh;
vhn    = vh;
Vn     = [uhn; vhn];
pn     = p;

% right hand side evaluations, initialized at zero
ku     = zeros(Nu,s_RK);
kv     = zeros(Nv,s_RK);
kp     = zeros(Np,s_RK);


% divergence operator
M      = options.discretization.M;
% boundary condition for divergence operator (known from last time step)
yMn    = options.discretization.yM;

if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn,options);
end

% 
% uLo_n  = uLo; uLe_n = uLe; uUp_n  = uUp; uRi_n = uRi;
% dudtLo_RK = zeros(Nx+1,s_RK); dudtUp_RK = zeros(Nx+1,s_RK);
% dudtLe_RK = zeros(Ny+1,s_RK); dudtRi_RK = zeros(Ny+1,s_RK);
% 
% uLo_i_n = uLo_i; uLe_i_n = uLe_i; uUp_i_n = uUp_i; uRi_i_n = uRi_i;
% dudtLo_RK_i = zeros(Nux_in,s_RK); dudtUp_RK_i = zeros(Nux_in,s_RK);
% dudtLe_RK_i = zeros(Ny,s_RK); dudtRi_RK_i = zeros(Ny,s_RK);
% 
% vLo_n  = vLo; vLe_n = vLe; vUp_n  = vUp; vRi_n = vRi;
% dvdtLo_RK = zeros(Nx+1,s_RK); dvdtUp_RK = zeros(Nx+1,s_RK);
% dvdtLe_RK = zeros(Ny+1,s_RK); dvdtRi_RK = zeros(Ny+1,s_RK);
% 
% vLo_i_n  = vLo_i; vLe_i_n = vLe_i; vUp_i_n = vUp_i; vRi_i_n = vRi_i;
% dvdtLo_RK_i = zeros(Nx,s_RK); dvdtUp_RK_i = zeros(Nx,s_RK);
% dvdtLe_RK_i = zeros(Nvy_in,s_RK); dvdtRi_RK_i = zeros(Nvy_in,s_RK);

ti = tn;

for i_RK=1:s_RK
    % at i=1 we calculate F_1, p_2 and u_2
    % ...
    % at i=s we calculate F_s, p_(n+1) and u_(n+1)
        
    
    % right-hand side for ti based on current velocity field uh, vh at
    % level i
    % this includes force evaluation at ti 
    % this exclu
    % boundary conditions should have been set through set_bc_vectors
    [~,Fu,Fv]  = F(uh,vh,p,ti,options);
    
    % store right-hand side of stage i
    ku(:,i_RK) = Fu + Omu_inv.*(Gx*p); % this removes the pressure contribution
    kv(:,i_RK) = Fv + Omv_inv.*(Gy*p); 
    
    % update velocity current stage by sum of F_i's until this stage,
    % weighted with Butcher tableau coefficients
    % this gives u_(i+1), and for i=s gives u_(n+1)
    utemp      = ku*A_RK(i_RK,:)';
    vtemp      = kv*A_RK(i_RK,:)';    
    Vtemp      = [utemp;vtemp];    
    
    % to make the velocity field u_(i+1) at t_(i+1) divergence-free we need 
    % the boundary conditions at t_(i+1)  
    ti         = tn + c_RK(i_RK)*dt;
    
    if (options.BC.BC_unsteady == 1)
        options = set_bc_vectors(ti,options);
    end
    yM   = options.discretization.yM;
%     y_px = options.discretization.y_px;
%     y_py = options.discretization.y_py;
    
    % divergence of R is directly calculated with M
    f       = (M*Vtemp + (yM-yMn)/dt)/c_RK(i_RK);
    
    % we should have sum(f) = 0 for periodic and no-slip BC
    % solve the Poisson equation for the pressure, but not for the first
    % step if the boundary conditions are steady
    if (options.BC.BC_unsteady==1 || i_RK>1)
        dp = pressure_poisson(f,ti,options);
    else
        dp = zeros(Np,1);
    end
%     p  = p + dp;    
    kp(:,i_RK) = dp;

   
    % update velocity current stage, which is now divergence free
    uh      = uhn + dt*(utemp - c_RK(i_RK)*Omu_inv.*(Gx*dp));
    vh      = vhn + dt*(vtemp - c_RK(i_RK)*Omv_inv.*(Gy*dp));

    
end

% gather solution into V
V = [uh;vh];


if (options.BC.BC_unsteady == 1)
    % for steady BC we skip this step and do an additional pressure solve 
    % that saves a pressure solve for i=1 in the next time step
    
    % make a suitable combination of pressures from stages that satisfy the
    % C(2) simplifying condition
    % three-stage method
%     p = -3*kp(:,2) +4*kp(:,3);
    % four-stage method
%     p = -2*kp(:,2) + 3*kp(:,4);

    % make a suitable combination with the theory of Hairer
%     W = inv(A_RK)*diag(c_RK);
%     p = kp*W(end,:)';

    % standard method
    p = kp(:,end);
else
    p = pressure_additional_solve(uh,vh,p,tn+dt,options);
end
