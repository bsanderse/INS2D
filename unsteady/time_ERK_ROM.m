%% general explicit Runge-Kutta method for ROM

% number of unknowns (modes) in ROM
M  = options.rom.M;

%% get coefficients of RK method 
% need to do this only once, as long as the RK method does not change in
% time
if (t==options.time.t_start)
    
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
    
    % vector with time instances
    c_RK = [c_RK(2:end);1]; % 1 is the time level of final step
    
end

%% preprocessing

% store variables at start of time step
tn     = t;
Rn     = R;
qn     = q;

% right hand side evaluations, initialized at zero
kR     = zeros(M,s_RK);
% array for the pressure
% kp     = zeros(Np,s_RK);
%

ti = tn;

for i_RK=1:s_RK
    % at i=1 we calculate F_1, p_2 and u_2
    % ...
    % at i=s we calculate F_s, p_(n+1) and u_(n+1)
    
    
    % right-hand side for ti based on current field R at
    % level i (this includes force evaluation at ti)
    % note that input q is not used in F_ROM if the basis is div-free
    [~,F_rhs]  = F_ROM(R,q,ti,options);
    
    % store right-hand side of stage i
    kR(:,i_RK) = F_rhs;
    
    % update coefficients R of current stage by sum of F_i's until this stage,
    % weighted with the Butcher tableau coefficients
    % this gives R_(i+1), and for i=s gives R_(n+1)
    Rtemp      = kR*A_RK(i_RK,:)';
    
    % time level of the computed stage
    ti         = tn + c_RK(i_RK)*dt;
    
    if (options.rom.div_free == 0)

        f  = (options.rom.Mdiv*(Rn/dt + Rtemp))/c_RK(i_RK);
        dq = pressure_poisson_ROM(f,t,options);

%         if (options.rom.precompute_pressure == 1)
            % approach 1: (with precomputed matrices)
%             Gq = options.rom.G*q;
%         elseif (options.rom.precompute_pressure == 0)
%             % approach 2: evaluate convection on FOM level, then map back
%             Gp = options.rom.B' * options.discretization.G * p;
%         end        
        
        R  = Rn + dt*(Rtemp - c_RK(i_RK)*(options.rom.G*dq));

    else
        % div-free ROM, no pressure needed
        % update ROM coefficients current stage
        R  = Rn + dt*Rtemp;
    end
    

    
end


if (options.rom.div_free == 1)
    if (options.rom.pressure_recovery == 1)
        q = pressure_additional_solve_ROM(R,tn+dt,options);
        p = getFOM_pressure(q,t,options);
    end
else
    % this might be improved like in FOM methods
    q = dq;
end