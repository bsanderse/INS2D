%% kinetic energy (actually mostly computations)
if options.verbosity.energy_verbosity == 1
    NV = options.grid.NV;
    Np = options.grid.Np;
    
    nu = 1/options.fluid.Re;
%     Q_h = options.discretization.Q_h;
%     
%     k_diff_ = - nu*norm(Q_h*V);
    
    D_h = options.discretization.D_h;
    k_diff_ = -nu*(V'*D_h*V);
    
    K_h = options.discretization.K_h;
    I_h = options.discretization.I_h;
    A_h = options.discretization.A_h;
    y_A = options.discretization.y_A;
    y_I = options.discretization.y_I;
    
    NF = numel(y_I);
    
    k_conv_ = -V'*(K_h*(spdiags(I_h*V+y_I,0,NF,NF)*(A_h*V+y_A)));
    
    p_h = zeros(Np,1);
    p_h = pressure_additional_solve(V,p_h,t,options);
    
    M_h = options.discretization.M;
    
    k_pres_ = (M_h*V)'*p_h;
    
    y_G = [options.discretization.y_px; ...
        options.discretization.y_py];
    
    y_D = [options.discretization.yDiffu; ...
        options.discretization.yDiffv];

    k_presBC_ = -V'*y_G;
    
    k_diffBC_ = nu*(V'*y_D);
    
    [Fx, Fy] = force(V,t,options,0);
    
    k_force_ = V'*[Fx;Fy];
    
    k_diff(n) = k_diff_;
    k_conv(n) = k_conv_;
    k_pres(n) = k_pres_;
    k_presBC(n) = k_presBC_;
    k_diffBC(n) = k_diffBC_;
    k_force(n) = k_force_;
    
%     profile off
%     profile viewer
%     17
end