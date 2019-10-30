% Runge-Kutta 4
    b_RK = [1;2;2;1]/6;
    a_RK = [1; 1; 2]/2;
    A_RK = spdiags(a_RK,-1,4,4);
    c_RK = [0; 1/2; 1/2; 1];
% Runge-Kutta 3
%     b_RK = [(8/15)-(17/60);0;3/4];
%     A_RK = spalloc(3,3,3);
%     A_RK(2,1) = 8/15;
%     A_RK(3,1) = (8/15)-(17/60);
%     A_RK(3,2) = 5/12;
%     c_RK = [0; A_RK(2,1); A_RK(3,1)+A_RK(3,2)];
% Runge-Kutta 2
%       b_RK = [1/2; 1/2];
%       A_RK = [0 0; 1 0];
%       c_RK = [0;1];
% Runge-Kutta 1 (Forward Euler)
%         b_RK = 1;
%         A_RK = 0;
%         c_RK = 0;
%     
% we work with the following 'shifted' Butcher tableau, because A_RK(1,1)
% is always zero for explicit methods
A_RK = [A_RK(2:end,:); b_RK'];

% number of stages
s_RK = length(c_RK);
c_RK = [c_RK;1]; % 1 is the time level of final step

% store variables at start of time step
tn     = t;
uhn    = uh;
vhn    = vh;
pn     = p;

ku     = zeros(Nu,s_RK);
kv     = zeros(Nv,s_RK);
yM1 = yM;
for i_RK=1:s_RK
    % at i=1 we calculate p_1, K_1, and u_2
    % ...
    % at i=s we calculate p_s, K_s, and u_(n+1)
    
    t       = tn + c_RK(i_RK)*dt;

    if (BC_unsteady == 1)
        boundary_conditions;
        interpolate_bc;
        operator_bc_momentum;
        operator_bc_divergence;
    end
    force;
    
    % convection
    cu      = uh;
    cv      = vh;
    convection;

    % diffusion
    d2u     = Diffu*uh + yDiffu;
    d2v     = Diffv*vh + yDiffv;
 
    % right hand side for pressure equation =
    % current conv+diffusion, without pressure
    % add pressure BC y_px and y_py?
    utemp   = Omu_inv.*( - du2dx - duvdy + d2u + Fx);
    vtemp   = Omv_inv.*( - duvdx - dv2dy + d2v + Fy);
    R       = [utemp;vtemp];
    % divergence of R is directly calculated with M
    f       = M*R + ydM;    
    
    % update yM
%     t       = tn + c_RK(i_RK+1)*dt;
%     boundary_conditions;
%     interpolate_bc;
%     operator_bc_divergence;
%    
%          
%     if(i_RK>1)
%     ftemp = [uhn + dt*ku(:,1:i_RK-1)*A_RK(i_RK,1:i_RK-1)' + dt*A_RK(i_RK,i_RK)*utemp;...
%              vhn + dt*kv(:,1:i_RK-1)*A_RK(i_RK,1:i_RK-1)' + dt*A_RK(i_RK,i_RK)*vtemp];
%     else
%     ftemp = [uhn + dt*A_RK(i_RK,i_RK)*utemp;...
%              vhn + dt*A_RK(i_RK,i_RK)*vtemp];
%     end
%     f      = 1/(dt*A_RK(i_RK,i_RK)) * (M*ftemp + yM);
  
    % we should have sum(f) = 0 for periodic and no-slip BC
    % solve the Poisson equation for the pressure                      %
    % using pre-determined LU decomposition
%     if (i_RK>1)
        pressure_poisson;
        p = dp;
%     end

% 
%     % update conv+diffusion with pressure gradient
    ku(:,i_RK) = utemp - Omu_inv.*(Gx*p + y_px);
    kv(:,i_RK) = vtemp - Omv_inv.*(Gy*p + y_py);
    
    % update velocity current stage
    uh      = uhn + dt*ku*A_RK(i_RK,:)';
    vh      = vhn + dt*kv*A_RK(i_RK,:)';

%     max(abs(M*[uh;vh]+yM))
%     keyboard;
end
V = [uh;vh];

t = tn; % time update is performed in solver_unsteady;

pressure_additional_solve;