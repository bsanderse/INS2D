function q = pressure_additional_solve_ROM(R,t,options)
% additional pressure solve
% returns ROM pressure coefficients; to get pressure use getROM_pressure(q)
        



    if (options.rom.pressure_precompute==0)        
        % without precomputing:
        Bp     = options.rom.Bp; 

        % note that time-varying BCs are not implemented yet
        V         = getFOM_velocity(R,t,options);
        % call F with additional argument nopressure=1 
        % to only get convection and diffusion
        % note that p is then irrelevant        
        [~,F_FOM] = F(V,V,0,t,options,0,1);
        Om_inv    = options.grid.Om_inv;
        f         = Bp' * options.discretization.M * (Om_inv.*F_FOM);
    
    elseif (options.rom.pressure_precompute == 1)
        % with precomputing (see operator_rom.m)        
%         if (options.case.force_unsteady == 1)
%             % for unsteady forcing, we have to compute
%             [Fx, Fy] = force(V,t,options,getJacobian);
%             options.rom.ppe_force = options.rom.P_PPE*[Fx;Fy];
%         end
        options.rom.ppe_force = options.rom.ppe_force*(1+sin(pi*t)); % see also F_ROM.m
        
        f   = options.rom.ppe_force + options.rom.ppe_bc + options.rom.ppe_linear*R + options.rom.ppe_quad*kron(R,R);
    
    end

    if (options.rom.pressure_mean == 1)
        % subtract mean pressure field
        f = f - options.rom.ppe_mean;
    end
        
        
    
    L   = options.rom.L;
    U   = options.rom.U;
    b   = L\f;
    % pressure coefficients of reduced order model    
    q   = U\b; 

end