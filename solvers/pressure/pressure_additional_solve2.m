function p = pressure_additional_solve2(V,t,disc_ydM,options)
% % additional pressure solve
% % make the pressure compatible with the velocity field. this should
% % also result in same order pressure as velocity

% compute pressure-like variable which ensures divergence-freeness

    Om_inv = options.grid.Om_inv;

%     % get updated BC for ydM
%     if (options.BC.BC_unsteady == 1)   
%         options = set_bc_vectors(t,options);
%     end
    
    M  = options.discretization.M;
    
    % note: F already contains G*p with the current p
    % we therefore effectively solve for the pressure difference
    [~,R,~] = F(V,V,0,t,options,0,1);
        
    f  = M*(Om_inv.*R) + disc_ydM;

%     Om = options.grid.Om;
%     Om_inv12 = Om.^-.5;
%     [Q1,~] = qr(Om_inv12.*M',0);
%     Q1t = Om_inv.^.5.*Q1;
%     f2 = M*Q1t*Q1t'*R + disc_ydM;
% 
%     dp2 = pressure_poisson(f2,t,options);
% norm(dp-dp2)
    
    dp = pressure_poisson(f,t,options);

    p  = dp;

end