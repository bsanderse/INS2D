% function [k_diff_,k_conv_,k_pres_,k_presBC_,k_diffBC_,k_force_] ...
function  [k_analysis,k_sum,p_h] = kinetic_energy_analysis(V,t,dt,options,k_analysis)
%% kinetic energy (actually mostly computations)
if options.verbosity.energy_verbosity == 1
    NV = options.grid.NV;
    Np = options.grid.Np;
    
    nu = 1/options.fluid.Re;
%     Q_h = options.discretization.Q_h;
%     
%     k_diff_ = - nu*norm(Q_h*V);
    
    D_h = options.discretization.D_h;
    k_diff_ = (V'*D_h*V);
    
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
    
    k_pres_ = +(M_h*V)'*p_h;
%     k_pres_ = 0; %+(M_h*V)'*p_h; %botch!!!
% warning('botch here!')
    
    y_G = [options.discretization.y_px; ...
        options.discretization.y_py];
    
    y_D = [options.discretization.yDiffu; ...
        options.discretization.yDiffv];

    k_presBC_ = -V'*y_G;
    
    k_diffBC_ = (V'*y_D);
    
    [Fx, Fy] = force(V,t,options,0);
    
    k_force_ = V'*[Fx;Fy];
    
    %%
    BC = options.BC;
    obc = max(strcmp({BC.u.left,BC.u.right,BC.u.low,BC.u.up},'mvp-obc'));
    if obc
        gO = @(V) options.BC.gO(V) + 0*V;
        gO_factor = options.grid.gO_factor;
        
        Conv_diag = options.grid.C;
        y_O = spdiags(V,0,NV,NV)*(Conv_diag*V)...
            - spdiags(gO_factor,0,NV,NV)*(spdiags(gO(V),0,NV,NV)*V);
        k_obc_ = V'*y_O;
    else
        k_obc_ = 0;
        y_O = 0;
    end
    %%
    
    k_analysis.k_diff = dt*k_diff_;
    k_analysis.k_conv = dt*k_conv_;
    k_analysis.k_pres = dt*k_pres_;
    k_analysis.k_presBC = dt*k_presBC_;
    k_analysis.k_diffBC = dt*k_diffBC_;
    k_analysis.k_force = dt*k_force_;
    k_analysis.k_obc = dt*k_obc_;
    
    k_sum = k_diff_ + k_conv_ + k_pres_ + k_force_ ...
          + k_diffBC_ + k_presBC_ + k_obc_;
    k_sum = 2*dt*k_sum;
    
%     F_rhs = - (K_h*(spdiags(I_h*V+y_I,0,NF,NF)*(A_h*V+y_A))) ...
%             + D_h*V + M_h'*p_h - y_G + y_D + [Fx;Fy] + y_O; 
        
%     F_rhs = - (K_h*(spdiags(I_h*V+y_I,0,NF,NF)*(A_h*V+y_A))) ...
%             + D_h*V - y_G + y_D + [Fx;Fy] + y_O;
%     B = (options.rom.B);
%     F_rhs1 = B'*F_rhs
% %     
%     R = getROM_velocity(V,t,options);
%     [~,F_rhs2] = F_ROM(R,0,t,options,0);
%     F_rhs2
%     norm(F_rhs1-F_rhs2)
%     17
    
%     k_analysis.k_diff(n) = k_diff_;
%     k_analysis.k_conv(n) = k_conv_;
%     k_analysis.k_pres(n) = k_pres_;
%     k_analysis.k_presBC(n) = k_presBC_;
%     k_analysis.k_diffBC(n) = k_diffBC_;
%     k_analysis.k_force(n) = k_force_;
    
%     profile off
%     profile viewer
%     17
end