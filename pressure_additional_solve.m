%% additional pressure solve
% make the pressure compatible with the velocity field. this should
% also result in same order pressure as velocity

if (p_add_solve==1)
    % evaluate BC and force at end of time step
    t = tn + dt;
    if (BC_unsteady == 1)   
        boundary_conditions;
        interpolate_bc;
        operator_bc_momentum;
        operator_bc_divergence;
    end
    force;
    t = tn;    
      
    % convection
    cu     = uh;
    cv     = vh;
    convection;

    % diffusion
    diffusion;

    Ru     =  - convu + d2u + Fx - y_px;    
    Rv     =  - convv + d2v + Fy - y_py;
    R      = Om_inv.*[Ru;Rv];
    f      = M*R + ydM;

    pressure_poisson;

    p = dp;
end