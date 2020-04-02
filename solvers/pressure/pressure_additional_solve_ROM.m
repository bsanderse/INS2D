function p = pressure_additional_solve_ROM(V,t,options)
% additional pressure solve
% make the pressure compatible with the velocity field. this should
% also result in same order pressure as velocity


    Bp     = options.rom.Bp; 
    Om_inv = options.grid.Om_inv;

    % call F with additional argument nopressure=1 
    % to only get convection and diffusion
    % note that p is then irrelevant
    
    % in the future, precompute stuff
    % note that time-varying BCs are not implemented yet
    [~,F_FOM] = F(V,V,0,t,options,0,1);
    f   = Bp' * options.discretization.M * (Om_inv.*F_FOM);
    L   = options.rom.L;
    U   = options.rom.U;
    b   = L\f;
    Q   = U\b; % coefficients of reduced order model
    p   = Bp*Q; % pressure at FOM level

end