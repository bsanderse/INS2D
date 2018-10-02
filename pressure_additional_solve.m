function p = pressure_additional_solve(uh,vh,p,t,options)
% additional pressure solve
% make the pressure compatible with the velocity field. this should
% also result in same order pressure as velocity

    Om_inv = options.grid.Om_inv;

    if (options.BC.BC_unsteady == 1)   
        options = set_bc_vectors(t,options);
    end
    
    M  = options.discretization.M;
    % note: derivative of BC!
    ydM = options.discretization.ydM;
    
    % note: F already contains G*p with the current p
    % we therefore effectively solve for the pressure difference
    [~,Ru,Rv] = F(uh,vh,p,t,options);
    
    R  = [Ru;Rv];
    
    f  = M*R + ydM;
    
    dp = pressure_poisson(f,t,options);

    p  = p+dp;

end