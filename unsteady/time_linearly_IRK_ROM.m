%% general implicit Runge-Kutta method for ROM

% number of unknowns (modes) in ROM
M  = options.rom.M;

%% get coefficients of RK method
if (t==options.time.t_start)
    
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
    
end

%% preprocessing

% store variables at start of time step
tn     = t;
Rn     = R;

% tj contains the time instances at all stages, tj = [t1;t2;...;ts]
tj    = tn + c_RK*dt;

% iteration counter
i = 0;
% iteration error
error_nonlinear = zeros(nonlinear_maxit,1);

% Vtot contains all stages and is ordered as [u1;v1;u2;v2;...;us;vs];
% initialize with the solution at tn
Rtotn = kron(ones(s_RK,1),Rn);

% index in global solution vector
indxR = 1:M*s_RK;

% starting guess for intermediate stages
Rj    = Rtotn;

Qj = Rj;

% initialize right-hand side for all stages
[~,F_rhs,~]  = F_multiple_ROM(Rj,[],tj,options,0);
% initialize momentum residual
fmom   = - (Rj - Rtotn)/dt + A_RK_ext*F_rhs;
% initialize residual
f      = fmom;

if (options.solversettings.nonlinear_Newton == 1) % approximate Newton
    % Jacobian based on current solution un 
    [~,~,Jn] = F_ROM(Rn,[],tn,options,1);
    % form iteration matrix, which is now fixed during iterations
    dfmom = (Om_sM/dt - kron(A_RK,Jn));
    %
    Z = dfmom;

end

while (max(abs(f))> options.solversettings.nonlinear_acc)
   
    if (options.solversettings.nonlinear_Newton == 1) 
        % approximate Newton   
        % do not rebuild Z
        dQj = Z\f;

    elseif (options.solversettings.nonlinear_Newton == 2)
        % full Newton
        [~,~,J] = F_multiple_ROM(Rj,[],tj,options,1);
        % form iteration matrix
        dfmom   = Om_sM/dt - A_RK_ext*J;
      
        Z  = dfmom;
           
        % get change
        dQj = Z\f;
    end
    
    % update solution vector
    Qj  = Qj + dQj;
    Rj  = Qj(indxR);
    
    % update iteration counter
    i   = i+1;    
    
    % evaluate rhs for next iteration and check residual based on
    % computed Rj
    [~,F_rhs,~] = F_multiple_ROM(Rj,[],tj,options,0);
    fmom   = - (Rj - Rtotn)/dt + A_RK_ext*F_rhs;
    
    f = fmom; 
    
    error_nonlinear(i) = max(abs(f));
    if (i>nonlinear_maxit)
        error(['Newton not converged in ' num2str(nonlinear_maxit) ' iterations']);
    end
    
end

nonlinear_its(n) = i;

% solution at new time step with b-coefficients of RK method
R = Rn + dt*(b_RK_ext*F_rhs);

if (options.rom.pressure_recovery == 1)
    q = pressure_additional_solve_ROM(R,tn+dt,options);
    p = getFOM_pressure(q,t,options);
end
