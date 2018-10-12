%% general implicit Runge-Kutta method

% (unsteady) Dirichlet boundary points are not part of solution vector but
% are prescribed in a 'strong' manner via the uBC and vBC functions

%% grid info
Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = Nu+Nv;
Np = options.grid.Np;
Om_inv = options.grid.Om_inv;
Om = options.grid.Om;


%% get coefficients of RK method
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
Om_sNV    = kron(Is,spdiags(Om,0,NV,NV));
A_RK_ext  = kron(A_RK,speye(NV));
b_RK_ext  = kron(b_RK',speye(NV));
c_RK_ext  = spdiags(c_RK,0,s_RK,s_RK);

%% preprocessing

% store variables at start of time step
tn     = t;
Vn     = V;
pn     = p;
qn     = [Vn; pn];

% tj contains the time instances at all stages, tj = [t1;t2;...;ts]
tj    = tn + c_RK*dt;

% gradient operator
G      = options.discretization.G;
Gtot   = kron(A_RK,G); % could also use 1 instead of c_RK and later scale the pressure
% divergence operator
M      = options.discretization.M;
Mtot   = kron(Is,M);
% finite volumes
Omtot  = kron(ones(s_RK,1),Om);

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

% zero block in iteration matrix
Z2 = spalloc(s_RK*Np,s_RK*Np,0);

% iteration counter
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
[~,F_rhs,~]  = F_multiple(Vj,pj,tj,options,0);
% initialize momentum residual
fmom   = - (Omtot.*Vj - Omtot.*Vtotn)/dt + A_RK_ext*F_rhs;
% initialize mass residual
fmass  = - (Mtot*Vj + yMtot);
f      = [fmom;fmass];

if (options.solversettings.nonlinear_Newton == 1) % approximate Newton
    % Jacobian based on current solution un 
    [~,~,Jn] = F(Vn,pn,tn,options,1);
    % form iteration matrix, which is now fixed during iterations
    dfmom = (Om_sNV/dt - kron(A_RK,Jn));
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
        dfmom   = Om_sNV/dt-A_RK_ext*J;
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
    [~,F_rhs,~] = F_multiple(Vj,pj,tj,options,0);
    fmom   = - (Omtot.*Vj - Omtot.*Vtotn)/dt + A_RK_ext*F_rhs;
    fmass  = - (Mtot*Vj + yMtot);
    
    f = [fmom; fmass];
    
    error_nonlinear(i) = max(abs(f));
    if (i>nonlinear_maxit)
        error(['Newton not converged in ' num2str(nonlinear_maxit) ' iterations']);
    end
    
end


% solution at new time step with b-coefficients of RK method
V = Vn + dt*Om_inv.*(b_RK_ext*F_rhs);

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
