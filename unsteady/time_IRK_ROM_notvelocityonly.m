function [Vnew,pnew,iterations,k_delta,k_analysis] = time_IRK_ROM_notvelocityonly(Rn,pn,tn,dt,options,k_analysis)

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

c_RK_ext_p  = kron(c_RK_ext,speye(Np));
c_RK_ext_V  = kron(c_RK_ext,speye(NV));
A_RK_ext_p  = kron(A_RK,speye(Np));

%% rom input
B = options.rom.B;
P = options.rom.P;
M_ROM = options.rom.M;

B_ext  = kron(speye(s_RK),B);
P_ext  = kron(speye(s_RK),P);

%% preprocessing

% store variables at start of time step
% tn     = t;
% Vn     = V;
% pn     = p;
% qn     = [Vn; pn];
% Rn = R;
% Vn = B*Rn;
% pn = pressure_additional_solve(V,p,t,options);

p  = pn;

% tj contains the time instances at all stages, tj = [t1;t2;...;ts]
tj    = tn + c_RK*dt;

% gradient operator
G      = options.discretization.G;
% Gtot   = kron(A_RK,G); % could also use 1 instead of c_RK and later scale the pressure
Gtot2  = kron(Is,G); % could also use 1 instead of c_RK and later scale the pressure
% divergence operator
M_h      = options.discretization.M;
Mtot   = kron(Is,M_h);
% finite volumes
Omtot  = kron(ones(s_RK,1),Om);
Om_invtot  = kron(ones(s_RK,1),Om_inv);
% Laplace operator
L = M_h*(Om_inv.*G);
Ltot = kron(Is,L);

% boundary condition for divergence operator
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn,options);
end
yMn    = options.discretization.yM;
yMn_tot = kron(Is,yMn);

% to make the velocity field u_(i+1) at t_(i+1) divergence-free we need 
% the boundary conditions at t_(i+1)  
if (options.BC.BC_unsteady == 1)
    yMtot = zeros(Np,s_RK);
    for i=1:s_RK
        ti      = tj(i);
        options = set_bc_vectors(ti,options);
        yMtot(:,i) = options.discretization.yM;            
    end
    yMtot = yMtot(:);
else
    yMtot = kron(ones(s_RK,1),yMn);
end

% zero block in iteration matrix
Z2 = spalloc(s_RK*Np,s_RK*Np,0);

% iteration counter
i = 0;
% iteration error
nonlinear_maxit = options.solversettings.nonlinear_maxit;
error_nonlinear = zeros(nonlinear_maxit,1);

% Vtot contains all stages and is ordered as [u1;v1;u2;v2;...;us;vs];
% initialize with the solution at tn
% Vtotn = kron(ones(s_RK,1),Vn);
Rtotn = kron(ones(s_RK,1),Rn);
ptotn = kron(ones(s_RK,1),pn);

% index in global solution vector
% indxV = 1:NV*s_RK;
% indxp = (NV*s_RK+1):(NV+Np)*s_RK;
indxR = 1:M_ROM*s_RK;
indxp = M_ROM*s_RK+1:(M_ROM+Np)*s_RK;

% starting guess for intermediate stages => this can be improved, see e.g.
% the Radau,Gauss4, or Lobatto scripts
% Vj    = Vtotn;
Rj    = Rtotn;
pj    = ptotn;
% Qj    = [Vj;pj];
% Rj    = P_ext*(Om.*Vj);
Qj    = [Rj;pj];

Vj = B_ext*Rj;
Vtotn = B_ext*Rn;

% % initialize right-hand side for all stages
% [~,F_rhs,~]  = F_multiple(Vj,Vj,pj,tj,options,0,1);
% % initialize momentum residual
% % fmom   = - (Omtot.*Vj - Omtot.*Vtotn)/dt + A_RK_ext*F_rhs;
% % initialize mass residual
% % fmass  = - (Mtot*Vj + yMtot);
% % Mass conservation does in general not hold for standard POD Galerkin ROM
% % Hence, use (7.1) (projected onto POD basis) and (7.2) in Sanderse PhD instead
% fmom   = - (Omtot.*Vj - Omtot.*Vtotn)/dt + A_RK_ext*F_rhs - c_RK_ext_V*Gtot2*pj;
% fmass = - c_RK_ext_p*Ltot*pj + A_RK_ext_p*Mtot*F_rhs - (yMtot-yMn_tot)/dt;
% % f      = [fmom;fmass];
% f      = [P_ext*fmom;fmass];

RHS_param.indxR = indxR;
RHS_param.indxp = indxp;
RHS_param.tj = tj;
RHS_param.options = options;
RHS_param.Omtot = Omtot;
RHS_param.Vtotn = Vtotn;
RHS_param.dt = dt;
RHS_param.A_RK_ext = A_RK_ext;
RHS_param.c_RK_ext_V = c_RK_ext_V;
RHS_param.Gtot2 = Gtot2;
RHS_param.c_RK_ext_p = c_RK_ext_p; 
RHS_param.Ltot = Ltot;
RHS_param.A_RK_ext_p = A_RK_ext_p;
RHS_param.Mtot = Mtot;
RHS_param.b = (yMtot-yMn_tot)/dt;
RHS_param.P_ext = P_ext;
RHS_param.B_ext = B_ext;

% norm(f-RHS(Qj,RHS_param))
    [f,F_rhs] = RHS(Qj,RHS_param);


if (options.solversettings.nonlinear_Newton == 1) % approximate Newton
    % Jacobian based on current solution un 
    [~,~,Jn] = F(Vn,Vn,pn,tn,options,1);
    % form iteration matrix, which is now fixed during iterations
%     dfmom = (Om_sNV/dt - kron(A_RK,Jn));
    dfmomdV = (Om_sNV/dt - kron(A_RK,Jn));
    dfmomdp = c_RK_ext_V*Gtot2;
    dfmasdV = -kron(A_RK,M_h*Jn);
    dfmasdp = c_RK_ext_p*Ltot;
    %
%     Z = [dfmom Gtot; ...
%          Mtot Z2];
%     Z = [P*dfmom*B P*Gtot; ...
%          Mtot*B Z2];
    Z = [P_ext*dfmomdV*B_ext P_ext*dfmomdp; dfmasdV*B_ext dfmasdp];
    % determine LU decomposition; often this is too slow
%     [L,U] = lu(Z);
end

while (max(abs(f))> options.solversettings.nonlinear_acc)
    
%     if options.verbosity.debug_mode == 1
%         max(abs(f))
%     end
    
    if (options.solversettings.nonlinear_Newton == 1) 
        % approximate Newton   
        % do not rebuild Z
        dQj = Z\f;
        % re-use the LU decomposition (often too slow):
%         dQj = U\(L\f);
        
        
    elseif (options.solversettings.nonlinear_Newton == 2)
        % full Newton
        [~,~,J] = F_multiple(Vj,Vj,pj,tj,options,1);
        % form iteration matrix
%         dfmom   = Om_sNV/dt-A_RK_ext*J;
%         Z   = [dfmom Gtot; ...
%                Mtot Z2];        
%         Z = [P*dfmom*B P*Gtot; ...
%             Mtot*B Z2];

    dfmomdV = Om_sNV/dt - A_RK_ext*J;
    dfmomdp = c_RK_ext_V*Gtot2;
    dfmasdV = -A_RK_ext_p*Mtot*J;
    dfmasdp = c_RK_ext_p*Ltot;
    %
%     Z = [dfmom Gtot; ...
%          Mtot Z2];
%     Z = [P*dfmom*B P*Gtot; ...
%          Mtot*B Z2];
    Z = [P_ext*dfmomdV*B_ext P_ext*dfmomdp; dfmasdV*B_ext dfmasdp];

        % get change
        dQj = Z\f;
    end
    
    % update solution vector
%     Qj  = Qj + dQj;
%     Vj  = Qj(indxV);
%     pj  = Qj(indxp);
    Qj  = Qj + dQj;
    Rj  = Qj(indxR);
    pj  = Qj(indxp);
    Vj  = B_ext*Rj;
    
    
    
    % update iteration counter
    i   = i+1;    
    
%     % evaluate rhs for next iteration and check residual based on
%     % computed Vj, pj
%     [~,F_rhs,~] = F_multiple(Vj,Vj,pj,tj,options,0,1);
% %     fmom   = - (Omtot.*Vj - Omtot.*Vtotn)/dt + A_RK_ext*F_rhs;
% %     fmass  = - (Mtot*Vj + yMtot);
% %     
% % %     f = [fmom; fmass];
% %     f      = [P*fmom;fmass];
%     
%     % Mass conservation does in general not hold for standard POD Galerkin ROM
%     % Hence, use (7.1) (projected onto POD basis) and (7.2) in Sanderse PhD instead
%     fmom   = - (Omtot.*Vj - Omtot.*Vtotn)/dt + A_RK_ext*F_rhs - c_RK_ext_V*Gtot2*pj;
%     fmass = - c_RK_ext_p*Ltot*pj + A_RK_ext_p*Mtot*F_rhs - (yMtot-yMn_tot)/dt;
%     % f      = [fmom;fmass];
%     f      = [P_ext*fmom;fmass];
% 
%     norm(f-RHS(Qj,RHS_param))
    [f,F_rhs] = RHS(Qj,RHS_param);

    error_nonlinear(i) = max(abs(f));
    if (i>nonlinear_maxit)
        Q1 = fsolve(@(Q) RHS(Q,RHS_param),Qj);
         error(['Newton not converged in ' num2str(nonlinear_maxit) ' iterations']);
    end
    
end

% store number of iterations
iterations = i;
iterations

% solution at new time step with b-coefficients of RK method
% R = Rn + dt*Om_inv.*(b_RK_ext*F_rhs);
R = Rn + P*(dt*Om_inv.*(b_RK_ext*F_rhs));

if options.verbosity.energy_verbosity == 1
%     k_delta = Vn'*(Om_inv.*(b_RK_ext*F_rhs))*dt;
%     k_delta = 2*dt*(Vn'*(b_RK_ext*F_rhs));
%     k_delta = 2*dt*(V'*(b_RK_ext*F_rhs));    
%     k_delta = 2*dt*(Vj'*(b_RK_ext*F_rhs));
    NV = options.grid.NV;
    Vjj = reshape(Vj,NV,numel(Vj)/NV);
    Fjj = reshape(F_rhs,NV,numel(Vj)/NV);
    k_delta = 2*dt*sum(b_RK'*(Vjj'*Fjj)); 
    
    % only for GL1 !!!
    [k_analysis,k_sum,p_h] = kinetic_energy_analysis(Vj,tj,dt,options,k_analysis);
%     k_delta - k_sum
%     norm(p_h-pj)
else
    k_delta = -666;
end

% make V satisfy the incompressibility constraint at n+1; this is only
% needed when the boundary conditions are time-dependent
% for stiffly accurate methods, this can also be skipped (e.g. Radau IIA) -
% this still needs to be implemented

if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn+dt,options);
    yM      = options.discretization.yM;
    V = B*R;
    f       = (1/dt)*(M_h*V + yM);

    dp      = pressure_poisson(f,tn+dt,options);

    V       = V - dt*Om_inv.*(G*dp);

end


if (options.BC.BC_unsteady == 1)
    
    if (options.solversettings.p_add_solve == 1)
        V = B*R;
        p = pressure_additional_solve(V,p,tn+dt,options);
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

Vnew = R;
pnew = p;


function [RHS_,F_rhs_] = RHS(Qj,RHS_param)
    indxR_ = RHS_param.indxR;
    indxp_ = RHS_param.indxp;
    tj_ = RHS_param.tj;
    options_ = RHS_param.options;
    Omtot_ = RHS_param.Omtot;
    Vtotn_ = RHS_param.Vtotn;
    dt_ = RHS_param.dt;
    A_RK_ext_ = RHS_param.A_RK_ext;
    c_RK_ext_V_ = RHS_param.c_RK_ext_V;
    Gtot2_ = RHS_param.Gtot2;
    c_RK_ext_p_ = RHS_param.c_RK_ext_p; 
    Ltot_ = RHS_param.Ltot;
    A_RK_ext_p_ = RHS_param.A_RK_ext_p;
    Mtot_ = RHS_param.Mtot;
    b_ = RHS_param.b;
    P_ext_ = RHS_param.P_ext;
    B_ext_ = RHS_param.B_ext;

    Rj_  = Qj(indxR_);
    pj_  = Qj(indxp_);

    Vj_ = B_ext_*Rj_;
    [~,F_rhs_,~]  = F_multiple(Vj_,Vj_,pj_,tj_,options_,0,1);

    % Mass conservation does in general not hold for standard POD Galerkin ROM
    % Hence, use (7.1) (projected onto POD basis) and (7.2) in Sanderse PhD instead
    fmom_   = - (Omtot_.*Vj_ - Omtot_.*Vtotn_)/dt_ ...
            + A_RK_ext_*F_rhs_ - c_RK_ext_V_*Gtot2_*pj_;
    fmass_ = - c_RK_ext_p_*Ltot_*pj_ + A_RK_ext_p_*Mtot_*F_rhs_ - b_;
    % f      = [fmom;fmass];
    f_      = [P_ext_*fmom_;fmass_];
    RHS_ = f_;
end

end