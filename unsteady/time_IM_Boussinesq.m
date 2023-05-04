function [Vnew,pnew,Tnew,iterations] = time_IM_Boussinesq(Vn,pn,Tn,tn,dt,options)

%% implicit midpoint method for the Boussinesq system

% (unsteady) Dirichlet boundary points are not part of solution vector but
% are prescribed in a 'strong' manner via the uBC and vBC functions

%% grid info
Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = Nu+Nv;
Np = options.grid.Np;
NVT = NV+Np; % number of unknowns velocity and temperature
Om_inv  = options.grid.Om_inv;
Om      = options.grid.Om;
Omp_inv = options.grid.Omp_inv;
Omp     = options.grid.Omp;
OmVT     = [Om;Omp];
OmVT_inv = [Om_inv;Omp_inv];
% mass matrix (with finite volume sizes)
Omtot = spdiags(OmVT,0,NVT,NVT);

%% get coefficients of RK method
% if (isnumeric(options.time.RK))
%     options.time.RK = num2str(options.time.RK);
% end
% implicit midpoint = gauss-legendre 1-stage
[A_RK,b_RK,c_RK] = getRKmethod('GL1');


%% preprocessing

% tj contains the time instance at the intermediate stage
tj    = tn + c_RK*dt;

% gradient operator
G      = options.discretization.G;
% divergence operator
M      = options.discretization.M;

% boundary condition for divergence operator at tn
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn,options);
end
yMn    = options.discretization.yM;

% to make the velocity field u_(i+1) at t_(i+1) divergence-free we need
% the boundary conditions at t_(i+1)
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tj,options);
    yMj   = options.discretization.yM;
else
    yMj = yMn;
end

% zero blocks in iteration matrix
Z2 = spalloc(Np,Np+Np,0);
Z3 = spalloc(Np,Np,0);
% extend G
Gext = [G; Z3];

% iteration counter
i = 0;
% iteration error
nonlinear_maxit = options.solversettings.nonlinear_maxit;
error_nonlinear = zeros(nonlinear_maxit,1);


% index in global solution vector; we order as [V; T; p];
indxV = (1:NV)';
indxT = (indxV(end)+1:indxV(end)+Np)';
indxp = (indxT(end)+1:indxT(end)+Np)';

% velocity and temperature unknowns
indxVT = [indxV;indxT];

% starting guess for intermediate stages => this can be improved, see e.g.
% the Radau, Gauss4, or Lobatto scripts
Vj    = Vn;
Tj    = Tn;
pj    = pn;

% total solution vector
Qn    = [Vn;Tn;pn];
Qj    = [Vj;Tj;pj];

% initialize right-hand side; this returns a vector
% [FV; FT]
[~,F_rhs,~]  = F(Vj,Vj,pj,Tj,tj,options,0);
% initialize mass residual
fmass  = - (M*Vj + yMj);
% initialize momentum and temperature residual
fmom_temp = - OmVT.*(Qj(indxVT) - Qn(indxVT))/dt + A_RK*F_rhs;

f    = [fmom_temp;fmass];


% we want to solve V1, T1, p1 from
% Omega*(V1-Vn)/dt = 0.5*FV(V1,T1) - 0.5*G*p1
% Omega*(T1-Tn)/dt = 0.5*FT(V1,T1)
% M*V1 + yM = 0
% and then update
% Omega*(Vn+1 - Vn)/dt = FV(V1,T1)
% Omega*(Tn+1 - Tn)/dt = FT(V1,T1)

if (options.solversettings.nonlinear_Newton == 1) % approximate Newton
    % Jacobian based on current solution un
    [~,~,Jn] = F(Vn,Vn,pn,Tn,tn,options,1);
    % form iteration matrix, which is now fixed during iterations
    dfmom_temp = (Omtot/dt - A_RK*Jn);
    %
    Z = [dfmom_temp A_RK*Gext; ...
         M Z2];
    % determine LU decomposition; often this is too slow
    %     [L,U] = lu(Z);
end

while (max(abs(f))> options.solversettings.nonlinear_acc)
    
    if (options.solversettings.nonlinear_Newton == 1)
        % approximate Newton
        % do not rebuild Z
        dQj = Z\f;
%         dQj = gmres(Z,f,[],1e-8,100);
        % re-use the LU decomposition (often too slow):
        %         dQj = U\(L\f);
        
        
    elseif (options.solversettings.nonlinear_Newton == 2)
        % full Newton
        [~,~,J] = F(Vj,Vj,pj,Tj,tj,options,1);
        % form iteration matrix
        dfmom_temp   = Omtot/dt - A_RK*J;
        Z   = [dfmom_temp A_RK*Gext; ...
            M Z2];
        % get change
        dQj = Z\f;
    end
    
    % update solution vector
    Qj  = Qj + dQj;
    Vj  = Qj(indxV);
    Tj  = Qj(indxT);
    pj  = Qj(indxp);
    
    % update iteration counter
    i   = i+1;
    
    % evaluate rhs for next iteration and check residual based on
    % computed Vj, pj
    [~,F_rhs,~] = F(Vj,Vj,pj,Tj,tj,options,0);
    fmom_temp   = - OmVT.*(Qj(indxVT) - Qn(indxVT))/dt + A_RK*F_rhs;
    fmass       = - (M*Vj + yMj);
    
    f = [fmom_temp; fmass];
    
    error_nonlinear(i) = max(abs(f));
    if (i>nonlinear_maxit)
        error(['Newton not converged in ' num2str(nonlinear_maxit) ' iterations']);
    end
    
end

% store number of iterations
iterations = i;
error_nonlinear(1:i);

% solution at new time step with b-coefficients of RK method
Q = Qn(indxVT) + dt*OmVT_inv.*(b_RK*F_rhs);
V = Q(indxV);
T = Q(indxT);

% make V satisfy the incompressibility constraint at n+1; this is only
% needed when the boundary conditions are time-dependent
% for stiffly accurate methods, this can also be skipped (e.g. Radau IIA) -
% this still needs to be implemented

if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn+dt,options);
    yM      = options.discretization.yM;
    f       = (1/dt)*(M*V + yM);
    
    dp      = pressure_poisson(f,tn+dt,options);
    
    V       = V - dt*Om_inv.*(G*dp);
    
end


if (options.BC.BC_unsteady == 1)
    
    if (options.solversettings.p_add_solve == 1)
        p = pressure_additional_solve(V,p,0,tn+dt,options);
    else
        
        % make a suitable combination of pressures from stages that satisfy the
        % C(2) simplifying condition
        % three-stage method
        %     p = -3*kp(:,2) +4*kp(:,3);
        % four-stage method
        %     p = -2*kp(:,2) + 3*kp(:,4);
        
        % make a suitable combination with the theory of Hairer
        %     W = inv(A_RK)*diag(c_RK);
        %     p = kp*W(end,:)';
        
        % standard method; take last pressure
        p = pj(end-Np+1:end);
    end
else
    % for steady BC we do an additional pressure solve
    % that saves a pressure solve for i=1 in the next time step
    %     p = pressure_additional_solve(V,p,tn+dt,options);
    
    % standard method; take pressure of last stage
    p = pj(end-Np+1:end);
    
end

Vnew = V;
pnew = p;
Tnew = T;