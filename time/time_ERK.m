% general explicit Runge-Kutta method
% integration of boundary values

Nu = options.grid.Nu;
Nv = options.grid.Nv;
Np = options.grid.Np;

Om_inv = options.grid.Om_inv;

G = options.discretization.G;


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

% number of stages
c_RK = [c_RK(2:end);1]; % 1 is the time level of final step

% store variables at start of time step
tn     = t;
uhn    = uh;
vhn    = vh;
Vn     = [uhn; vhn];
pn     = p;

% right hand side evaluations, initialized at zero
k     = zeros(Nu+Nv,s_RK);
% k     = zeros(Nv,s_RK);
kp     = zeros(Np,s_RK);


% divergence operator
M      = options.discretization.M;
% boundary condition for divergence operator (known from last time step)
yMn    = options.discretization.yM;
yM     = yMn;    


if (options.BC.BC_unsteady == 1)
    options = set_bc_vectors(tn,options);
end

ti = tn;

for i_RK=1:s_RK
    % at i=1 we calculate F_1, p_2 and u_2
    % ...
    % at i=s we calculate F_s, p_(n+1) and u_(n+1)
        
    
    % right-hand side for ti based on current velocity field uh, vh at
    % level i
    % this includes force evaluation at ti 
    % boundary conditions should have been set through set_bc_vectors
    [~,F_rhs]  = F(uh,vh,p,ti,options);
    
    % store right-hand side of stage i
    % we remove the pressure contribution Gx*p and Gy*p (but not the
    % vectors y_px and y_py)
    k(:,i_RK) = F_rhs + Om_inv.*(G*p);
    
    % update velocity current stage by sum of F_i's until this stage,
    % weighted with Butcher tableau coefficients
    % this gives u_(i+1), and for i=s gives u_(n+1)
    Vtemp      = k*A_RK(i_RK,:)';
    
    % to make the velocity field u_(i+1) at t_(i+1) divergence-free we need 
    % the boundary conditions at t_(i+1)  
    ti         = tn + c_RK(i_RK)*dt;
    if (options.BC.BC_unsteady == 1)
        options = set_bc_vectors(ti,options);
        yM      = options.discretization.yM;    
    end
    
    % divergence of intermediate velocity field is directly calculated with M
    f       = (M*Vtemp + (yM-yMn)/dt)/c_RK(i_RK);
    
    % we should have sum(f) = 0 for periodic and no-slip BC
    % solve the Poisson equation for the pressure, but not for the first
    % step if the boundary conditions are steady
    if (options.BC.BC_unsteady==1 || i_RK>1)
        % the time ti below is only for output writing
        dp = pressure_poisson(f,ti,options);
    else
        dp = pn;
    end
    % store pressure
    kp(:,i_RK) = dp;
   
    % update velocity current stage, which is now divergence free
    V  = Vn + dt*(Vtemp - c_RK(i_RK)*Om_inv.*(G*dp));
    uh = V(1:Nu);
    vh = V(Nu+1:Nu+Nv);
    
end

% gather solution into V
V = [uh;vh];


if (options.BC.BC_unsteady == 1)
    % for steady BC we skip this step and do an additional pressure solve 
    % that saves a pressure solve for i=1 in the next time step
    
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
else
    p = pressure_additional_solve(uh,vh,p,tn+dt,options);
end
