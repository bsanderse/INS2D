function [Vnew,pnew] = time_ERK(Vn,pn,tn,dt,options)

%% general explicit Runge-Kutta method

% (unsteady) Dirichlet boundary points are not part of solution vector but
% are prescribed in a 'strong' manner via the uBC and vBC functions

%% grid info
Nu = options.grid.Nu;
Nv = options.grid.Nv;
Np = options.grid.Np;
Om_inv = options.grid.Om_inv;

%% get coefficients of RK method
% make character string if necessary
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

%% preprocessing
% store variables at start of time step
% tn     = t;
% Vn     = V;
% pn     = p;
V = Vn;
p = pn;

% right hand side evaluations, initialized at zero
kV     = zeros(Nu+Nv,s_RK);
% array for the pressure
kp     = zeros(Np,s_RK);

if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn,options);
end

% gradient operator
G      = options.discretization.G;
% divergence operator
M      = options.discretization.M;
% boundary condition for divergence operator (known from last time step)
yMn    = options.discretization.yM;
yM     = yMn;

ti = tn;

%% start looping over stages

for i_RK=1:s_RK
    % at i=1 we calculate F_1, p_2 and u_2
    % ...
    % at i=s we calculate F_s, p_(n+1) and u_(n+1)
    
    
    % right-hand side for ti based on current velocity field uh, vh at
    % level i
    % this includes force evaluation at ti and pressure gradient
    % boundary conditions will be set through set_bc_vectors inside F
    % the pressure p is not important here, it will be removed again in the
    % next step
    [~,F_rhs]  = F(V,V,p,ti,options);
    
    % store right-hand side of stage i
    % by adding G*p we effectively REMOVE the pressure contribution Gx*p and Gy*p (but not the
    % vectors y_px and y_py)
    kV(:,i_RK)  = Om_inv.*(F_rhs + G*p);
    
    % update velocity current stage by sum of F_i's until this stage,
    % weighted with Butcher tableau coefficients
    % this gives u_(i+1), and for i=s gives u_(n+1)
    Vtemp      = kV*A_RK(i_RK,:)';
    
    % to make the velocity field u_(i+1) at t_(i+1) divergence-free we need
    % the boundary conditions at t_(i+1)
    ti         = tn + c_RK(i_RK)*dt;
    if (options.BC.BC_unsteady == 1)
        options = set_bc_vectors(ti,options);
        yM      = options.discretization.yM;
    end
    
    % divergence of intermediate velocity field is directly calculated with M
    % old formulation:
%     f       = (M*Vtemp + (yM-yMn)/dt)/c_RK(i_RK);
    % new formulation, prevents growth of constraint errors:
    % instead of -yMn we use +M*Vn; they are the same up to machine
    % precision but using the latter prevents error accumulation
    f       = (M*(Vn/dt+Vtemp) + yM/dt)/c_RK(i_RK);
    % note: we should have sum(f) = 0 for periodic and no-slip BC
    
    % solve the Poisson equation for the pressure, but not for the first
    % step if the boundary conditions are steady
    if (options.BC.BC_unsteady==1 || i_RK>1)
        % the time ti below is only for output writing
        dp = pressure_poisson(f,ti,options);
    else % BC steady AND i_RK=1
        dp = pn;
    end
    % store pressure
    kp(:,i_RK) = dp;
    
    % update velocity current stage, which is now divergence free
    V  = Vn + dt*(Vtemp - c_RK(i_RK)*Om_inv.*(G*dp));
    
end


if (options.BC.BC_unsteady == 1)
    
    if (options.solversettings.p_add_solve == 1)
        p = pressure_additional_solve(V,p,tn+dt,options);
    else
        
        % make a suitable combination of pressures from stages that satisfy the
        % C(2) simplifying condition
        % three-stage method
        %     p = -3*kp(:,2) +4*kp(:,3);
        % four-stage method
        %     p = -2*kp(:,2) + 3*kp(:,4);
        
        % make a suitable combination with the theory of Hairer
        %     W = inv(A_RK)*diag(c_RK);
        %     p = kp*W(end,:)';
        
        % standard method
        p = kp(:,end);
    end
else
    % for steady BC we do an additional pressure solve
    % that saves a pressure solve for i=1 in the next time step
    p = pressure_additional_solve(V,p,tn+dt,options);
end

Vnew = V;
pnew = p;
