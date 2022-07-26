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

% right hand side evaluations, initialized at zero
kR     = zeros(M,s_RK);
if options.rom.bc_recon == 2
    NV = options.grid.NV;
    kR_FOM  = zeros(NV,s_RK);
end
% array for the pressure
% kp     = zeros(Np,s_RK);
%

ti = tn;

options = set_bc_vectors(t,options);
yMn      = options.discretization.yM;

for i_RK=1:s_RK
    % at i=1 we calculate F_1, p_2 and u_2
    % ...
    % at i=s we calculate F_s, p_(n+1) and u_(n+1)
    
    if options.rom.bc_recon == 2
        V = getFOM_velocity(R,t,options);
        [~,F_rhs_FOM] = F(V,V,p,t,options,0,1);

        P = options.rom.P;
        F_rhs = P*F_rhs_FOM;
    elseif (options.rom.bc_recon == 4)
        [~,F_rhs]  = F_ROM_notvelocityonly(R,p,ti,options);
    else
        % right-hand side for ti based on current field R at
        % level i (this includes force evaluation at ti)
        % note that input p is not used in F_ROM
        [~,F_rhs]  = F_ROM(R,p,ti,options);
    end


    % store right-hand side of stage i
    kR(:,i_RK) = F_rhs;

    % update coefficients R of current stage by sum of F_i's until this stage,
    % weighted with the Butcher tableau coefficients
    % this gives R_(i+1), and for i=s gives R_(n+1)
    Rtemp      = kR*A_RK(i_RK,:)';
    if options.rom.bc_recon == 2
        %     if false
        kR_FOM(:,i_RK) = F_rhs_FOM;
        Rtemp_FOM = kR_FOM*A_RK(i_RK,:)';
    end

    % time level of the computed stage
    ti         = tn + c_RK(i_RK)*dt;
    if options.rom.bc_recon == 2

        if (options.BC.BC_unsteady == 1) || i_RK==1
            options = set_bc_vectors(ti,options);
            yM      = options.discretization.yM;
            ydM     = options.discretization.ydM;
            phi_bc  = options.rom.phi_bc;
        end
        % divergence of intermediate velocity field is directly calculated with M
        % old formulation:
        %     f       = (M*Vtemp + (yM-yMn)/dt)/c_RK(i_RK);
        % new formulation, prevents growth of constraint errors:
        % instead of -yMn we use +M*Vn; they are the same up to machine
        % precision but using the latter prevents error accumulation

        B = options.rom.B;
        % divergence operator
        M_h      = options.discretization.M;

        %         f       = (M_h*B*(Rn/dt+Rtemp) + yM/dt)/c_RK(i_RK);
        %         f       = (M_h*B*Rtemp + (yM-yMn)/dt)/c_RK(i_RK);
        if options.rom.time_discB4pres_elim == 1
            yM_diff = (yM-yMn)/dt;
        else
            yM_diff = ydM;
        end
        f       = (M_h*(Om_inv.*Rtemp_FOM) + yM_diff)/c_RK(i_RK);

        %         f       = (M_h*B*R + (yM-yMn)/dt)/c_RK(i_RK); %botch
        %         f       = (M_h*V + (yM-yMn)/dt)/c_RK(i_RK); %botch
        % note: we should have sum(f) = 0 for periodic and no-slip BC
        if options.verbosity.debug_mode==1
            norm(yM-yMn)
        end

        % solve the Poisson equation for the pressure, but not for the first
        % step if the boundary conditions are steady
        if (options.BC.BC_unsteady==1 || i_RK>1)
            % the time ti below is only for output writing
            dp = pressure_poisson(f,ti,options);
        else % BC steady AND i_RK=1
            %             dp = pn;
            dp = pressure_poisson(f,ti,options); % not efficient but safe
        end
        %         dp1 = dp;
        %         dp = pressure_additional_solve(B*Rtemp,p,ti,options);


        if options.verbosity.debug_mode == 1
            vis_p(dp,options)
        end
        % store pressure
        kp(:,i_RK) = dp;

        % gradient operator
        G_h      = options.discretization.G;

        % update ROM coefficients current stage
        P = options.rom.P;
        R  = Rn + dt*(Rtemp - P*(c_RK(i_RK)*Om_inv.*(G_h*dp)));
        %         R  = Rn + dt*(P*Rtemp - P*(c_RK(i_RK)*Om_inv.*(G_h*dp)));
    else

        % update ROM coefficients current stage
        R  = Rn + dt*Rtemp;
    end

    if options.rom.bc_recon == 5
        hatM = options.rom.hatM;
%         yMn = options.discretization.yM; % actually should use tilde yM
        % correct:
        phi_bc = options.rom.phi_bc;

%         yBCn = get_bc_vector_yBC(tn,options); % again wrong
        yBCn = phi_bc*get_a_bc(tn,options);
        yMn = get_yM(options,yBCn);

%         options = set_bc_vectors(ti,options); % actually should use tilde yM
%         yMi = options.discretization.yM; % actually should use tilde yM
        % correct:
%         yBCi = get_bc_vector_yBC(ti,options); % again wrong
        yBCi = phi_bc*get_a_bc(ti,options);
        yMi = get_yM(options,yBCi);


        Bp = options.rom.Bp;
%         norm(M_h*B*Rn+yMn)
%         norm(hatM*Rn + Bp'*yMn)

% norm(yMn-yMi)

        f = (hatM*(Rn/dt + Rtemp) + Bp'*yMi/dt)/c_RK(i_RK); % yM has opposite sign compared to theory
%         f = (hatM*R/dt + Bp'*yMi/dt)/c_RK(i_RK);

        hatG = - hatM';
        hatL = hatM*hatG;
        options.rom.hatL = hatL;
        q = hatL\f;

        R = R - dt*c_RK(i_RK)*(hatG*q); % = R - hatG*inv(hatL)*(hatM R+hat yM)
%         norm(hatM*R + Bp'*yMn) % = hatM R - hatM R - hat yM + hat yM = 0
%         norm(hatM*R + Bp'*yMi)

%% pressure hunt
% f1 = (hatM*Rtemp + Bp'*(yMi - yMn)/dt)/c_RK(i_RK);
% norm(f1-f)
%  M_h = options.discretization.M;
%  Om = options.grid.Om;
% 
%  [Q_,R_] = qr(Om.^-.5.*(M_h'),0);
%  Qt = Om.^-.5.*Q;
% 
%  f2 = Qt'*(R_temp+(yMi-yMn)/dt)/c_RK(i_RK);
%  L2 = Qt'*M_h
  

%%
    end

end


if (options.rom.pressure_recovery == 1)
    q = pressure_additional_solve_ROM(R,tn+dt,options);
    p = getFOM_pressure(q,t,options);
end