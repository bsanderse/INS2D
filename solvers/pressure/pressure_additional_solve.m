function p = pressure_additional_solve(V,p,t,options)
% additional pressure solve
% make the pressure compatible with the velocity field. this should
% also result in same order pressure as velocity

    Om_inv = options.grid.Om_inv;

    % get updated BC for ydM
    if (options.BC.BC_unsteady == 1)   
        options = set_bc_vectors(t,options);
    end
    
    M  = options.discretization.M;
    % note: time derivative of BC in ydM
    ydM = options.discretization.ydM;
    
    % note: F already contains G*p with the current p
    % we therefore effectively solve for the pressure difference
    [~,R,~] = F(V,V,p,t,options);
        
    f  = M*(Om_inv.*R) + ydM;
    
    dp = pressure_poisson(f,t,options);

    p  = p + dp;

end