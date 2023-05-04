function [Vnew,pnew,Tnew,iterations] = time_IM_Boussinesq_split(Vn,pn,Tn,tn,dt,options)
% conv_old are the convection terms of t^(n-1)
% output includes convection terms at t^(n), which will be used in next time step in
% the Adams-Bashforth part of the method

%% implicit midpoint method for the Boussinesq system,
% solved by iterating at each time step


%% grid info

Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = Nu+Nv;
Np = options.grid.Np;

indu = options.grid.indu;
indv = options.grid.indv;

Omu_inv = options.grid.Omu_inv;
Omv_inv = options.grid.Omv_inv;
Om_inv  = options.grid.Om_inv;
Omp_inv = options.grid.Omp_inv;
OmT     = options.grid.Omp;
OmV     = options.grid.Om;

% index in global solution vector; we order as [V; T; p];
indxV = (1:NV)';
indxT = (indxV(end)+1:indxV(end)+Np)';
% indxp = (indxT(end)+1:indxT(end)+Np)';

% velocity and temperature unknowns
% indxVT = [indxV;indxT];

%% coefficients of the method
% get coefficients of RK method
% if (isnumeric(options.time.RK))
%     options.time.RK = num2str(options.time.RK);
% end
% implicit midpoint = gauss-legendre 1-stage
[A_RK,b_RK,c_RK] = getRKmethod('GL1');

% Adams-Bashforth coefficients
% alfa1 = 3/2;
% alfa2 = -1/2;

% CN coefficients: should match with A_RK!
% theta = 1/2;
theta  = options.time.theta;
if (theta~=A_RK)
    error('theta value should be consistent with implicit midpoint value');
end


% tj contains the time instance at the intermediate stage
tj    = tn + c_RK*dt;


%% preprocessing
% store variables at start of time step
% tn     = t;
% Vn     = V;
% pn     = p;

% gradient operator
G      = options.discretization.G;
% divergence operator
M      = options.discretization.M;

% uh     = Vn(indu);
% vh     = Vn(indv);



%% convection from previous time step
% convu_old = rhs_terms_old.conv(indu);
% convv_old = rhs_terms_old.conv(indv);

%% evaluate BC and force at middle point
% assume for now that force is not depending on V
% unsteady BC at current time
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tj,options);
end

yMj   = options.discretization.yM;

% diffusion 
Diffu  = options.discretization.Diffu;
Diffv  = options.discretization.Diffv;
yDiffu = options.discretization.yDiffu;
yDiffv = options.discretization.yDiffv;

Gx   = options.discretization.Gx;
Gy   = options.discretization.Gy;
y_px = options.discretization.y_px;
y_py = options.discretization.y_py;

% diffusion of temperature 
switch options.case.boussinesq
    case 'temp'
        yDiffT = options.discretization.yDiffT;
end



%% starting guess for intermediate stages => this can be improved, see e.g.
% the Radau, Gauss4, or Lobatto scripts
Vj    = Vn;
Tj    = Tn;
pj    = pn;
dpj   = 0*pj;

% total solution vector
% Qn    = [Vn;Tn;pn];
% Qj    = [Vj;Tj;pj];

% initialize momentum and temperature residual
[~,F_rhs,~] = F(Vj,Vj,pj,Tj,tj,options,0);
ftemp   = - OmT.*(Tj - Tn)/dt + A_RK*F_rhs(indxT);
fmom    = - OmV.*(Vj - Vn)/dt + A_RK*F_rhs(indxV);
fmass   = - (M*Vj + yMj);

f    = [ftemp;fmom;fmass];

% iteration counter
i = 0;
% iteration error
nonlinear_maxit = options.solversettings.nonlinear_maxit;
error_nonlinear = zeros(nonlinear_maxit,1);


while (max(abs(f))> options.solversettings.nonlinear_acc)
    
    
    %% temperature equation 
    % first solve temperature equation with previous velocity fields
    switch options.case.boussinesq
        
        
        case 'temp'
            
            %
%             NT      = options.grid.NT;
%             DiffT   = options.discretization.DiffT;
            Omp_inv = options.grid.Omp_inv;
            
            % right-hand side of temperature equation
            convT   = convection_temperature(Tj,Vj,tj,options,0);
%             convT_old = rhs_terms_old.convT;
            
            FT    = Tn + A_RK*dt*Omp_inv.*( -convT + yDiffT );
            
            switch options.temp.incl_dissipation
                case 1
                    %  add dissipation to internal energy equation
                    Phi = dissipation(Vj,tj,options,0);
                    % the computed dissipation is basically V'*D*V, which has
                    % alfa1 as scaling
                    % in the internal energy equation we need alfa3, so we
                    % divide by gamma
%                     Phi_old = rhs_terms_old.Phi;
                    gamma  = options.temp.gamma;
                    FT  = FT + A_RK*dt*Omp_inv.*((1/gamma)*Phi);
            end
            
            % matrix arising from implicit diffusion: I-dt*Om_inv*D
            % Tnew = (speye(NT) - theta*dt*spdiags(Omp_inv,0,NT,NT)*DiffT) \ FT;
            % solve system for new temperature
            %  using pre-determined decomposition
            if (verLessThan('matlab','9.3'))
                b    = options.discretization.L_diffT\FT;
                Tj   = options.discretization.U_diffT\b;
            else
                Tj   = options.discretization.DiffT_impl\FT;
            end
            
            % add effect of temperature into momentum equation
            Fv_Temp  = options.discretization.AT_v*Tj;
            
        otherwise
            Tj = 0;
            
    end
    
    %% solve momentum equation
    % convection of current solution
    [convu, convv] = convection(Vj,Vj,tj,options,0);
    % forcing
    [Force_x,Force_y] = force(Vj,tj,options,0);
    switch options.case.boussinesq
        case 'temp'
            Force_y = Force_y + Fv_Temp;
    end
    
    % right hand side of the momentum equation update
    Rur = Vn(indu) + A_RK*Omu_inv*dt.*(- convu + ...
        + yDiffu + ...
        + Force_x - Gx*pj - y_px);
    
    Rvr = Vn(indv) + A_RK*Omv_inv*dt.*(- convv  + ...
        + yDiffv + ...
        + Force_y - Gy*pj - y_py);    

    
    % LU decomposition of diffusion part has been calculated already in
    % operator_convection_diffusion
    
    if (options.solversettings.poisson_diffusion==1)
        b   = options.discretization.L_diffu\Rur;
        Ru  = options.discretization.U_diffu\b;
        
        b   = options.discretization.L_diffv\Rvr;
        Rv  = options.discretization.U_diffv\b;
        
    elseif (options.solversettings.poisson_diffusion==3)
        [Ru,iter,norm1,norm2]=cg(L_diffu,int64(dia_diffu),int64(ndia_diffu),Rur,...
            CG_acc,int64(Nu),Ru,int64(CG_maxit));
        %                         iter
        [Rv,iter,norm1,norm2]=cg(L_diffv,int64(dia_diffv),int64(ndia_diffv),Rvr,...
            CG_acc,int64(Nv),Rv,int64(CG_maxit));
        %                         iter
    end
    
    Vtemp = [Ru; Rv];
    
        
    % boundary condition for the difference in pressure between time
    % steps; only non-zero in case of fluctuating outlet pressure
    y_dp = zeros(Nu+Nv,1);
    
    %% approach 1: use inverse of (I-0.5*dt*Om_inv*D) in computing the
    % generalized Poisson matrix; this works fine, but requires
    % precomputation of the generalized Poisson matrix, which is expensive: see
    % operator_convection_diffusion.m
% 
%     f    = (M*Vtemp + yMj)/(A_RK*dt) - M*y_dp; 
%     dpj  = options.solver_settings.A_split\f;
%     b    = options.discretization.L_diffu\(Omu_inv.*(Gx*dpj));
%     Ru_p = options.discretization.U_diffu\b;        
%     b    = options.discretization.L_diffv\(Omv_inv.*(Gy*dpj));
%     Rv_p = options.discretization.U_diffv\b;  
%     Vj   = Vtemp - A_RK*dt*[Ru_p;Rv_p];
% 
%     pj   = pj + dpj;
         
    % end approach 1
    
%  approach 2: iterate with Poisson matrix, put generalized Poisson matrix right hand side 
    A    = options.discretization.A;
    b    = options.discretization.L_diffu\(Omu_inv.*(Gx*dpj));
    Ru_p = options.discretization.U_diffu\b;
        
    b    = options.discretization.L_diffv\(Omv_inv.*(Gy*dpj));
    Rv_p = options.discretization.U_diffv\b;
    
    f    = (M*Vtemp + yMj)/(A_RK*dt) - M*y_dp - (M*[Ru_p;Rv_p] - A*dpj);
    
    % solve the Poisson equation for the pressure
    dpj  = pressure_poisson(f,tj,options);
    pj   = pj + dpj;

    
    % update velocity field
    b     = options.discretization.L_diffu\(Omu_inv.*(Gx*dpj));
    Ru_p2 = options.discretization.U_diffu\b;
        
    b     = options.discretization.L_diffv\(Omv_inv.*(Gy*dpj));
    Rv_p2 = options.discretization.U_diffv\b;
    
    Vj   = Vtemp - A_RK*dt*[Ru_p2;Rv_p2];

% end approach 2    

    
    % update iteration counter
    i   = i+1;
    
    % evaluate rhs and check residual based on
    % computed Vj, pj
    [~,F_rhs,~] = F(Vj,Vj,pj,Tj,tj,options,0);
    ftemp   = - OmT.*(Tj - Tn)/dt + A_RK*F_rhs(indxT);
    fmom    = - OmV.*(Vj - Vn)/dt + A_RK*F_rhs(indxV);
    fmass   = - (M*Vj + yMj);
    
    f = [ftemp; fmom; fmass];
    error_nonlinear(i) = max(abs(f));
    if (i>nonlinear_maxit)
        error(['Newton not converged in ' num2str(nonlinear_maxit) ' iterations']);
    end
    
end

% store number of iterations
iterations = i;
error_nonlinear(1:i);

% solution at new time step with b-coefficients of RK method
T = Tn + dt*Omp_inv.*(b_RK*F_rhs(indxT));
V = Vn + dt*Om_inv.*(b_RK*F_rhs(indxV));

% make V satisfy the incompressibility constraint at n+1; this is only
% needed when the boundary conditions are time-dependent
% for stiffly accurate methods, this can also be skipped (e.g. Radau IIA) -
% this still needs to be implemented
if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn+dt,options);
    yM      = options.discretization.yM;
    f       = (1/dt)*(M*V + yM);
    
    dp      = pressure_poisson(f,tn+dt,options);
    
    V       = V - dt*Om_inv.*(G*dp);
    
end
 
if (options.BC.BC_unsteady == 1)
    
    if (options.solversettings.p_add_solve == 1)
        p = pressure_additional_solve(V,pj,0,tn+dt,options);
    else       
        % standard method; take last pressure
        p = pj;
    end
else  
    % standard method; take pressure of last stage
    p = pj;    
end


% output:
Vnew = V;
pnew = p;
Tnew = T;

% % output convection at t^(n), to be used in next time step
% rhs_terms.conv = [convu;convv];
% 
% % get the temperature terms that are treated explicit in time
% switch options.case.boussinesq
%     
%     case 'temp'
%         
%         rhs_terms.convT = convT;
%         
%         switch options.temp.incl_dissipation
%             
%             case 1
%                 
%                 rhs_terms.Phi = Phi;
%         end
% end

