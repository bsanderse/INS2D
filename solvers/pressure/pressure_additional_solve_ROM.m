function q = pressure_additional_solve_ROM(R,t,options)
% additional pressure solve
% input: R, ROM coefficients of velocity field
% returns ROM pressure coefficients; to get pressure use getROM_pressure(q)


if (options.rom.pressure_precompute==0)
    % without precomputing:
    Bp     = options.rom.Bp;
    
    % note that time-varying BCs are not implemented yet
    V         = getFOM_velocity(R,t,options);
    % call F with additional argument nopressure=1
    % to only get convection and diffusion
    % note that p is then irrelevant
    [~,F_FOM] = F(V,V,0,0,t,options,0,1);
    Om_inv    = options.grid.Om_inv;
    f         = Bp' * options.discretization.M * (Om_inv.*F_FOM);
    
elseif (options.rom.pressure_precompute == 1)
    % with precomputing, see operator_rom.m for projection matrices
    
    if (options.force.force_unsteady == 1)        
        if (strcmp(options.case.project,'actuator_ROM'))
            if( options.time.t_start == t)
                % give warning only once
                warning('scaling unsteady force - actuator disk test case only!');
            end
            options.rom.ppe_force = options.rom.ppe_force*(1+sin(pi*t)); % see also F_ROM.m
        else
            error('you have unsteady forcing with precomputation: check if settings are correct');
        end
    end
    
    f   = options.rom.ppe_force + options.rom.ppe_bc + options.rom.ppe_linear*R + options.rom.ppe_quad*kron(R,R);
    
end

if (options.rom.pressure_mean == 1)
    % subtract mean pressure field
    f = f - options.rom.ppe_mean;
end

% pressure coefficients of reduced order model
%     % using LU decomposition:
%     L   = options.rom.L;
%     U   = options.rom.U;
%     b   = L\f;
%     q   = U\b;
q = options.rom.A_decomp \ f;

end