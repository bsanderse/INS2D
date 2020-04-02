%% general implicit Runge-Kutta method

% (unsteady) Dirichlet boundary points are not part of solution vector but
% are prescribed in a 'strong' manner via the uBC and vBC functions

% number of unknowns (modes) in ROM
M  = options.rom.M;
% 
% Nu = options.grid.Nu;
% Nv = options.grid.Nv;
% NV = Nu+Nv;
% Np = options.grid.Np;
% Om_inv = options.grid.Om_inv;
% Om = options.grid.Om;


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
Om_sM     = kron(Is,spdiags(ones(M,1),0,M,M));
A_RK_ext  = kron(A_RK,speye(M));
b_RK_ext  = kron(b_RK',speye(M));
c_RK_ext  = spdiags(c_RK,0,s_RK,s_RK);

%% preprocessing

% store variables at start of time step
tn     = t;
Rn     = R;
% Vn     = V;
% pn     = p;
% qn     = [Vn; pn];

% tj contains the time instances at all stages, tj = [t1;t2;...;ts]
tj    = tn + c_RK*dt;

% gradient operator
% G      = options.discretization.G;
% Gtot   = kron(A_RK,G); % could also use 1 instead of c_RK and later scale the pressure
% % divergence operator
% M      = options.discretization.M;
% Mtot   = kron(Is,M);
% finite volumes
% Omtot  = kron(ones(s_RK,1),ones(M,1)));
% 
% % boundary condition for divergence operator
% if (options.BC.BC_unsteady == 1)
%     options = set_bc_vectors(tn,options);
% end
% yMn    = options.discretization.yM;
% 
% % to make the velocity field u_(i+1) at t_(i+1) divergence-free we need 
% % the boundary conditions at t_(i+1)  
% if (options.BC.BC_unsteady == 1)
%     yMtot = zeros(Np,s_RK);
%     for i=1:s_RK
%         ti      = tj(i);
%         options = set_bc_vectors(ti,options);
%         yMtot(:,i) = options.discretization.yM;            
%     end
%     yMtot = yMtot(:);
% else
%     yMtot = kron(ones(s_RK,1),yMn);
% end

% zero block in iteration matrix
% Z2 = spalloc(s_RK*Np,s_RK*Np,0);

% iteration counter
i = 0;
% iteration error
error_nonlinear = zeros(nonlinear_maxit,1);

% Vtot contains all stages and is ordered as [u1;v1;u2;v2;...;us;vs];
% initialize with the solution at tn
Rtotn = kron(ones(s_RK,1),Rn);
% ptotn = kron(ones(s_RK,1),pn);

% index in global solution vector
indxR = 1:M*s_RK;
% indxp = (NV*s_RK+1):(NV+Np)*s_RK;

% starting guess for intermediate stages
Rj    = Rtotn;
% pj    = ptotn;
% Qj    = [Vj;pj];
Qj = Rj;

% initialize right-hand side for all stages
[~,F_rhs,~]  = F_multiple_ROM(Rj,[],tj,options,0);
% initialize momentum residual
fmom   = - (Rj - Rtotn)/dt + A_RK_ext*F_rhs;
% initialize mass residual
% fmass  = - (Mtot*Vj + yMtot);
% f      = [fmom;fmass];
f      = fmom;

if (options.solversettings.nonlinear_Newton == 1) % approximate Newton
    % Jacobian based on current solution un 
    [~,~,Jn] = F_ROM(Rn,[],tn,options,1);
    % form iteration matrix, which is now fixed during iterations
    dfmom = (Om_sM/dt - kron(A_RK,Jn));
    %
%     Z = [dfmom Gtot; ...
%          Mtot Z2];
    Z = dfmom;
    % determine LU decomposition; often this is too slow
%     [L,U] = lu(Z);
end

while (max(abs(f))> options.solversettings.nonlinear_acc)
   
    if (options.solversettings.nonlinear_Newton == 1) 
        % approximate Newton   
        % do not rebuild Z
        dQj = Z\f;
        % re-use the LU decomposition (often too slow):
%         dQj = U\(L\f);
        
        
    elseif (options.solversettings.nonlinear_Newton == 2)
        % full Newton
        [~,~,J] = F_multiple_ROM(Rj,[],tj,options,1);
        % form iteration matrix
        dfmom   = Om_sM/dt - A_RK_ext*J;
%         Z   = [dfmom Gtot; ...
%                Mtot Z2];        
        Z  = dfmom;
           
        % get change
        dQj = Z\f;
    end
    
    % update solution vector
    Qj  = Qj + dQj;
    Rj  = Qj(indxR);
%     pj  = Qj(indxp);
    
    % update iteration counter
    i   = i+1;    
    
    % evaluate rhs for next iteration and check residual based on
    % computed Vj, pj
    [~,F_rhs,~] = F_multiple_ROM(Rj,[],tj,options,0);
    fmom   = - (Rj - Rtotn)/dt + A_RK_ext*F_rhs;
%     fmass  = - (Mtot*Vj + yMtot);
    
    f = fmom; %; fmass];
    
    error_nonlinear(i) = max(abs(f));
    if (i>nonlinear_maxit)
        error(['Newton not converged in ' num2str(nonlinear_maxit) ' iterations']);
    end
    
end

nonlinear_its(n) = i;


% solution at new time step with b-coefficients of RK method
R = Rn + dt*(b_RK_ext*F_rhs);

% map back from reduced space to full order model space
V = B*R + Vbc;

if (options.rom.pressure_recovery == 1)
    p = pressure_additional_solve_ROM(V,t,options);
end
