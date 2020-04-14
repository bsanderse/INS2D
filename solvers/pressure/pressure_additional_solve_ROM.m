function q = pressure_additional_solve_ROM(R,t,options)
% additional pressure solve
% make the pressure compatible with the velocity field. this should
% also result in same order pressure as velocity


    Bp     = options.rom.Bp; 

    % call F with additional argument nopressure=1 
    % to only get convection and diffusion
    % note that p is then irrelevant
        
    if (options.rom.pressure_precompute==0)        
        % without precomputing:
        % note that time-varying BCs are not implemented yet
        V         = getFOM_velocity(R,t,options);
        [~,F_FOM] = F(V,V,0,t,options,0,1);
        Om_inv    = options.grid.Om_inv;
        f         = Bp' * options.discretization.M * (Om_inv.*F_FOM);
    
    elseif (options.rom.pressure_precompute == 1)
        % with precomputing:
        [~,F_rhs]  = F_ROM(R,0,t,options);
        f   = options.rom.M_ROM * F_rhs;
    
    end

    L   = options.rom.L;
    U   = options.rom.U;
    b   = L\f;
    % pressure coefficients of reduced order model    
    q   = U\b; 

end