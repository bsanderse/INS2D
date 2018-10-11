%% general implicit Runge-Kutta method

% (unsteady) Dirichlet boundary points are not part of solution vector but
% are prescribed in a 'strong' manner via the uBC and vBC functions

% grid info
Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = Nu+Nv;
Np = options.grid.Np;
Om_inv = options.grid.Om_inv;

% get coefficients of RK method
if (isnumeric(options.time.RK))
    options.time.RK = num2str(options.time.RK);
end
[A_RK,b_RK,c_RK] = getRKmethod(options.time.RK);
% RK_order = check_orderconditions(A_RK,b_RK,c_RK);
% number of stages
s_RK = length(b_RK);

options.time.A_RK = A_RK;
options.time.b_RK = b_RK;
options.time.c_RK = c_RK;
options.time.s_RK = s_RK;

% extend the Butcher tableau 
Is        = speye(s_RK);
I_sNV     = kron(Is,speye(NV));
A_RK_ext  = kron(A_RK,speye(NV));
b_RK_ext  = kron(b_RK',speye(NV));
c_RK_ext  = spdiags(c_RK,0,s_RK,s_RK);

% store variables at start of time step
tn     = t;
Vn     = V;
pn     = p;
qn     = [Vn; pn];

% tj contains the time instances at all stages, tj = [t1;t2;...;ts]
tj    = tn + c_RK*dt;

% right hand side evaluations, initialized at zero
% k      = zeros(NV,s_RK);
% array for the pressure
% kp     = zeros(Np,s_RK);




% gradient operator
G      = options.discretization.G;
Gtot   = kron(A_RK,spdiags(Om_inv,0,NV,NV)*G); % could also use 1 instead of c_RK and later scale the pressure
% divergence operator
M      = options.discretization.M;
Mtot   = kron(Is,M);

% boundary condition for divergence operator
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn,options);
end
yMn    = options.discretization.yM;

% to make the velocity field u_(i+1) at t_(i+1) divergence-free we need 
% the boundary conditions at t_(i+1)  
if (options.BC.BC_unsteady == 1)
    for i=1:s_RK
        ti      = tj(i);
        options = set_bc_vectors(ti,options);
        yMtot   = options.discretization.yM;            
    end
else
    yMtot = kron(ones(s_RK,1),yMn);
end

% zero block
Z2 = spalloc(s_RK*Np,s_RK*Np,0);

%
i = 0;

% capital U contains all stages and is ordered as [u1;v1;w1;...;u2;v2;w2;...;us;vs;ws];
% initialize U with the solution at tn
Vtotn    = kron(ones(s_RK,1),Vn);
ptotn    = kron(ones(s_RK,1),pn);
% Qn    = kron(ones(s_RK,1),qn);

% index in global solution vector
indxV = 1:NV*s_RK;
indxp = (NV*s_RK+1):(NV+Np)*s_RK;


% starting guess for intermediate stages
Vj    = Vtotn;
pj    = ptotn;
Qj    = [Vj;pj];

% initialize right-hand side for all stages
[maxres,F_rhs,~]  = F_multiple(Vj,pj,tj,options,0);
% initialize momentum residual
fmom   = - (Vj - Vtotn)/dt + A_RK_ext*F_rhs;
% initialize mass residual
fmass  = - (Mtot*Vj + yMtot);
f = [fmom;fmass];

if (options.solversettings.nonlinear_Newton == 1) % approximate Newton
    % Jacobian based on current solution un 
    [~,~,Jn] = F(Vn,pn,tn,options,1);
    % form iteration matrix, which is now fixed during iterations
    dfmom = (I_sNV/dt - kron(A_RK,Jn));
    %
    Z = [dfmom Gtot; ...
         Mtot Z2];
    % determine LU decomposition
    [L,U] = lu(Z);
end

while (max(abs(f))> options.solversettings.nonlinear_acc)
   
    if (options.solversettings.nonlinear_Newton == 1) 
        % approximate Newton        
        % re-use the LU decomposition
        dQj = U\(L\f);
        
    elseif (options.solversettings.nonlinear_Newton == 2)
        % full Newton
        [~,~,J] = F_multiple(Vj,pj,tj,options,1);
        % form iteration matrix
        dfmom   = I_sNV/dt-A_RK_ext*J;
        Z   = [dfmom Gtot; ...
               Mtot Z2];        
        % get change
        dQj = Z\f;
    end
    
    % update solution vector
    Qj  = Qj + dQj;
    Vj  = Qj(indxV);
    pj  = Qj(indxp);
    
    % update iteration counter
    i   = i+1;    
    
    % evaluate rhs for next iteration and check residual based on
    % computed Vj, pj
    [maxres,F_rhs,~] = F_multiple(Vj,pj,tj,options,0);
    fmom   = - (Vj - Vtotn)/dt + A_RK_ext*F_rhs;
    fmass  = - (Mtot*Vj + yMtot);
    
    f = [fmom; fmass];
    
    error_nonlinear(i) = max(abs(f));
%     max(abs(res));
    if (i>nonlinear_maxit)
        error(['Newton not converged in ' num2str(nonlinear_maxit) ' iterations']);
    end
    
    % right-hand side for ti based on current velocity field uh, vh at
    % level i
    % this includes force evaluation at ti 
    % boundary conditions should have been set through set_bc_vectors
%     [~,F_rhs]  = F_rhs(V,p,ti,options);
    
    % store right-hand side of stage i
    % we remove the pressure contribution Gx*p and Gy*p (but not the
    % vectors y_px and y_py)
%     k(:,i_RK)  = F_rhs + Om_inv.*(G*p);
    
    % update velocity current stage by sum of F_i's until this stage,
    % weighted with Butcher tableau coefficients
    % this gives u_(i+1), and for i=s gives u_(n+1)
%     Vtemp      = k*A_RK(i_RK,:)';
    

    
    % divergence of intermediate velocity field is directly calculated with M
%     f       = (M*Vtemp + (yM-yMn)/dt)/c_RK(i_RK);
%     
%     % we should have sum(f) = 0 for periodic and no-slip BC
%     % solve the Poisson equation for the pressure, but not for the first
%     % step if the boundary conditions are steady
%     if (options.BC.BC_unsteady==1 || i_RK>1)
%         % the time ti below is only for output writing
%         dp = pressure_poisson(f,ti,options);
%     else
%         dp = pn;
%     end
%     % store pressure
%     kp(:,i_RK) = dp;
%    
%     % update velocity current stage, which is now divergence free
%     V  = Vn + dt*(Vtemp - c_RK(i_RK)*Om_inv.*(G*dp));

end


% solution at new time step with b-coefficients of RK method
V = Vn + dt*b_RK_ext*F_rhs;

if (options.BC.BC_unsteady == 1)
    % make V satisfy the incompressibility constraint at n+1
%     t = tn + dt;

    options = set_bc_vectors(tn+dt,options);
    yM      = options.discretization.yM;
    f       = (1/dt)*(M*V + yM);

    dp      = pressure_poisson(f,tn+dt,options);

    V       = V - dt*Om_inv.*(G*dp);
    uh      = V(1:Nu);
    vh      = V(Nu+1:Nu+Nv);

end


if (options.BC.BC_unsteady == 0)
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
%     p = kp(:,end);
%     p = dp;
else
    p = pressure_additional_solve(V,p,tn+dt,options);
end
