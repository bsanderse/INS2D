% Runge-Kutta 4
%     b_RK = [1;2;2;1]/6;
%     a_RK = [1; 1; 2]/2;
%     A_RK = spdiags(a_RK,-1,4,4);
%     c_RK = [0; 1/2; 1/2; 1];
% RK4 satisfying C(2) for i=3
% A_RK=[0 0 0 0; 1/4 0 0 0; 0 1/2 0 0; 1 -2 2 0];
% b_RK=[1/6;0;2/3;1/6];
% c_RK=[0;1/4;1/2;1];
% RK4 satisfying C(2) for i=3 and c2=c3
A_RK=[0 0 0 0; 1/2 0 0 0; 1/4 1/4 0 0; 0 -1 2 0];
b_RK=[1/6;0;2/3;1/6];
c_RK=[0;1/2;1/2;1];
% RK4
% A_RK = [0 0 0 0; 1/3 0 0 0; -1/3 1 0 0; 1 -1 1 0];
% b_RK = [1/8;3/8;3/8;1/8];
% c_RK = [0;1/3;2/3;1];
% RK4 satisfying the second order condition for the pressure (but not third
% order)
% A_RK = [0 0 0 0; 1 0 0 0; 3/8 1/8 0 0; -1/8 -3/8 3/2 0];
% b_RK = [1/6; -1/18; 2/3; 2/9];
% c_RK = [0;1;1/2;1];
% Runge-Kutta 3 (Wray)
%     b_RK = [(8/15)-(17/60);0;3/4];
%     A_RK = zeros(3,3);
%     A_RK(2,1) = 8/15;
%     A_RK(3,1) = (8/15)-(17/60);
%     A_RK(3,2) = 5/12;
%     c_RK = [0; A_RK(2,1); A_RK(3,1)+A_RK(3,2)];
% RK3 satisfying C(2) for i=3
% A_RK = [0 0 0; 2/3 0 0; 1/3 1/3 0];
% b_RK = [1/4; 0; 3/4];
% c_RK = [0; 2/3; 2/3];
% RK3 satisfying the second order condition for the pressure
% A_RK = [0 0 0; 1/3 0 0; -1 2 0];
% b_RK = [0; 3/4; 1/4];
% c_RK = [0; 1/3; 1];
% Runge-Kutta 2
%       b_RK = [1/2; 1/2];
%       A_RK = [0 0; 1 0];
%       c_RK = [0;1];
% Runge-Kutta 1 (Forward Euler)
%         b_RK = 1;
%         A_RK = 0;
%         c_RK = 0;
%     
% Runge-Kutta 5th order, 6 stage
% A_RK = [0 0 0 0 0 0; 1/4 0 0 0 0 0; 1/8 1/8 0 0 0 0; 0 0 1/2 0 0 0; ...
%         3/16 -3/8 3/8 9/16 0 0; -3/7 8/7 6/7 -12/7 8/7 0];
% b_RK = [7/90;0;16/45;2/15;16/45;7/90];
% c_RK = [0;1/4;1/4;1/2;3/4;1];

% we work with the following 'shifted' Butcher tableau, because A_RK(1,1)
% is always zero for explicit methods
A_RK = [A_RK(2:end,:); b_RK'];

% number of stages
s_RK = length(c_RK);
c_RK = [c_RK(2:end);1]; % 1 is the time level of final step

% store variables at start of time step
tn     = t;
uhn    = uh;
vhn    = vh;
Vn     = [uhn; vhn];
pn     = p;

ku     = zeros(Nu,s_RK);
kv     = zeros(Nv,s_RK);
kp     = zeros(Np,s_RK);

yMn    = yM;

for i_RK=1:s_RK
    % at i=1 we calculate F_1, p_2 and u_2
    % ...
    % at i=s we calculate F_s, p_(n+1) and u_(n+1)
    
    % convection of U(i-1)
    cu      = uh;
    cv      = vh;
    convection;

    % diffusion of U(i-1)
    diffusion;
    % sum of F_i's until this stage
    % add pressure BC y_px and y_py??
    ku(:,i_RK) = Omu_inv.*( - convu + d2u + Fx);
    kv(:,i_RK) = Omv_inv.*( - convv + d2v + Fy);
    % update velocity current stage
    utemp      = ku*A_RK(i_RK,:)';
    vtemp      = kv*A_RK(i_RK,:)';    


    R       = [utemp;vtemp];
    
    t       = tn + c_RK(i_RK)*dt;
    
    if (i_RK==1 || c_RK(i_RK)~=c_RK(i_RK-1))
        % Matlab first checks the first condition, and if it is true it
        % skips the second condition, so i_RK-1 is not evaluated for i_RK=1
        % only necessary if c-coefficient is different from previous
        % c-coefficient    
        if (BC_unsteady == 1)

            boundary_conditions;
            interpolate_bc;     
            operator_bc_divergence;
            operator_bc_momentum;        
        end
        
        force;
        
    end
    
    % divergence of R is directly calculated with M
    f       = (1/c_RK(i_RK))*(M*R + (yM-yMn)/dt);
    
    % we should have sum(f) = 0 for periodic and no-slip BC
    % solve the Poisson equation for the pressure, but not for the first
    % step if the boundary conditions are steady
    if (BC_unsteady==1 || i_RK>1)
        pressure_poisson;
    else % BC steady AND i_RK=1
        dp = pn;
    end
    % store pressure    
    kp(:,i_RK) = dp;
   
    % update velocity current stage, which is now divergence free
    uh      = uhn + dt*(utemp - c_RK(i_RK)*Omu_inv.*(Gx*dp + y_px));
    vh      = vhn + dt*(vtemp - c_RK(i_RK)*Omv_inv.*(Gy*dp + y_py));

end
V = [uh;vh];

% for steady BC we skip this step and do an additional pressure solve 
% that saves a pressure solve for i=1 in the next time step
    
if (BC_unsteady == 1)

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
    p  = kp(:,end);
else
    p_add_solve = 1;
%     p = kp(:,end);
end

t = tn; % time update is performed in solver_unsteady;

pressure_additional_solve;