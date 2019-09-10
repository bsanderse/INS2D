%% general explicit Runge-Kutta method

% (unsteady) Dirichlet boundary points are not part of solution vector but
% are prescribed in a 'strong' manner via the uBC and vBC functions

% number of unknowns (modes) in ROM
M  = options.rom.M;
% Nu = options.grid.Nu;
% Nv = options.grid.Nv;
% Np = options.grid.Np;
% Om_inv = options.grid.Om_inv;

% get coefficients of RK method
if (isnumeric(options.time.RK))
    options.time.RK = num2str(options.time.RK);
end
[A_RK,b_RK,c_RK] = getRKmethod(options.time.RK);
% RK_order = check_orderconditions(A_RK,b_RK,c_RK);

% number of stages
s_RK = length(b_RK);

% we work with the following 'shifted' Butcher tableau, because A_RK(1,1)
% is always zero for explicit methods
A_RK = [A_RK(2:end,:); b_RK'];

% vector with time instances 
c_RK = [c_RK(2:end);1]; % 1 is the time level of final step

% store variables at start of time step
tn     = t;
Rn     = R;
% Vn     = V;
% pn     = p;

% right hand side evaluations, initialized at zero
kR     = zeros(M,s_RK);
% array for the pressure
% kp     = zeros(Np,s_RK);
% 
% if (options.BC.BC_unsteady == 1)
%     options = set_bc_vectors(tn,options);
% end

% gradient operator
% G      = options.discretization.G;
% divergence operator
% M      = options.discretization.M;
% boundary condition for divergence operator (known from last time step)
% yMn    = options.discretization.yM;
% yM     = yMn;    

ti = tn;

for i_RK=1:s_RK
    % at i=1 we calculate F_1, p_2 and u_2
    % ...
    % at i=s we calculate F_s, p_(n+1) and u_(n+1)
        
    
    % right-hand side for ti based on current field R at
    % level i (this includes force evaluation at ti)
    % boundary conditions will be set through set_bc_vectors inside F
    [~,F_rhs]  = F_ROM(R,p,ti,options);
    
    % store right-hand side of stage i
    kR(:,i_RK)  = F_rhs; 
    
    % update velocity current stage by sum of F_i's until this stage,
    % weighted with Butcher tableau coefficients
    % this gives u_(i+1), and for i=s gives u_(n+1)
    Rtemp      = kR*A_RK(i_RK,:)';
    
    % to make the velocity field u_(i+1) at t_(i+1) divergence-free we need 
    % the boundary conditions at t_(i+1)  
    ti         = tn + c_RK(i_RK)*dt;
%     if (options.BC.BC_unsteady == 1)
%         options = set_bc_vectors(ti,options);
%         yM      = options.discretization.yM;    
%     end
%     
%     % divergence of intermediate velocity field is directly calculated with M
%     f       = (M*Vtemp + (yM-yMn)/dt)/c_RK(i_RK);
%     
%     % we should have sum(f) = 0 for periodic and no-slip BC
%     % solve the Poisson equation for the pressure, but not for the first
%     % step if the boundary conditions are steady
%     if (options.BC.BC_unsteady==1 || i_RK>1)
%         % the time ti below is only for output writing
%         dp = pressure_poisson(f,ti,options);
%     else
%         dp = pn;
%     end
%     % store pressure
%     kp(:,i_RK) = dp;
   
    % update velocity current stage, which is now divergence free
    R  = Rn + dt*Rtemp; % - c_RK(i_RK)*Om_inv.*(G*dp));

end

% map back from reduced space to full order model space
V = B*R + Vbc;
