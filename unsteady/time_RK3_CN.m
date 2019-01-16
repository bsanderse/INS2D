% Runge-Kutta 3 a la Wray (see e.g. paper Knikker, IJNMF 2009)

% coefficients
beta1  = 4/15; beta2  =  1/15; beta3  = 1/6;
gamma1 = 0;    gamma2 =-17/60; gamma3 = -5/12;
sigma1 = 8/15; sigma2 =  5/12; sigma3 = 3/4;

dtn    = dt;
tn     = t;
uhn    = uh;
vhn    = vh;
pn     = p;

A_RKc  = [0 0 0 0; sigma1 0 0 0; sigma1+gamma2 sigma2 0 0];
b_RKc  = [sigma1+gamma2; sigma2+gamma3; sigma3; 0];

A_RKd  = [0 0 0 0; beta1 beta1 0 0; beta1 beta1+beta2 beta2 0];
b_RKd  = [beta1; beta1+beta2; beta2+beta3; beta3];

% A_RKd  = [0 0];
% b_RKd  = [1/2; 1/2];
% 
% A_RKc  = [0 0];
% b_RKc  = [1;0];

% assume the method has a_si = b_i and a_1i = 0
A_RKc  = [A_RKc(2:end,:); b_RKc'];
A_RKd  = [A_RKd(2:end,:); b_RKd'];

c_RK   = sum(A_RKc,2);

% number of stages
s_RK   = length(c_RK);

kuc    = zeros(Nu,s_RK); kud     = zeros(Nu,s_RK);
kvc    = zeros(Nv,s_RK); kvd     = zeros(Nv,s_RK);
kp     = zeros(Np,s_RK);

Z1     = spalloc(Np,Np,0);


for i_RK=1:s_RK
    
    % at i=1 we calculate F_1, p_2 and u_2
    % ...
    % at i=s we calculate F_s, p_(n+1) and u_(n+1)
    
    % convection of U(i-1)
    cu      = uh;
    cv      = vh;
    convection;

    % diffusion of U(i-1)
    d2u     = Diffu*uh + yDiffu;
    d2v     = Diffv*vh + yDiffv;

    % sum of F_i's until this stage
    kuc(:,i_RK) =  - du2dx - duvdy;
    kvc(:,i_RK) =  - duvdx - dv2dy;

    kud(:,i_RK) = d2u + Fx + y_px;
    kvd(:,i_RK) = d2v + Fy + y_py;
    
    t          = tn + c_RK(i_RK)*dt;
    
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
    
    % update velocity current stage, explicit part; 
    % use initial pressure to get second order
    utemp     = Omu.*uhn + dt*(kuc*A_RKc(i_RK,1:s_RK)' + kud*A_RKd(i_RK,1:s_RK)' + A_RKd(i_RK,i_RK+1)*(yDiffu + Fx));% - c_RK(i_RK)*dt*Gx*pn;
    vtemp     = Omv.*vhn + dt*(kvc*A_RKc(i_RK,1:s_RK)' + kvd*A_RKd(i_RK,1:s_RK)' + A_RKd(i_RK,i_RK+1)*(yDiffv + Fy));% - c_RK(i_RK)*dt*Gy*pn;     
    
%     % decoupled approach:
    utemp     = (spdiags(Omu,0,Nu,Nu) - A_RKd(i_RK,i_RK+1)*dt*Diffu)\utemp;
    vtemp     = (spdiags(Omv,0,Nv,Nv) - A_RKd(i_RK,i_RK+1)*dt*Diffv)\vtemp;
    R         = [utemp;vtemp];
    
%     divergence of R is directly calculated with M
    f         = 1/(dt*c_RK(i_RK))*(M*R + yM);
    
%     we should have sum(f) = 0 for periodic and no-slip BC
%     solve the Poisson equation for the pressure
    pressure_poisson;
    kp(:,i_RK) = dp;
   
%     update velocity current stage, which is now divergence free
    uh        = utemp - dt*c_RK(i_RK)*Omu_inv.*(Gx*dp); % + y_px
    vh        = vtemp - dt*c_RK(i_RK)*Omv_inv.*(Gy*dp); % + y_py

    p          = pn + dp;
    kp(:,i_RK) = dp;
    

    %     solve pressure and velocity coupled; more expensive but also more
    %     accurate
%     utemp     = Omu.*uhn + dt*(kuc*A_RKc(i_RK,1:s_RK)' + kud*A_RKd(i_RK,1:s_RK)' + A_RKd(i_RK,i_RK+1)*(yDiffu + Fx));
%     vtemp     = Omv.*vhn + dt*(kvc*A_RKc(i_RK,1:s_RK)' + kvd*A_RKd(i_RK,1:s_RK)' + A_RKd(i_RK,i_RK+1)*(yDiffv + Fy));    
%     
%     CD = [spdiags(Omu,0,Nu,Nu)-A_RKd(i_RK,i_RK+1)*dt*Diffu spalloc(Nu,Nv,0); ...
%           spalloc(Nv,Nu,0) spdiags(Omv,0,Nv,Nv)-A_RKd(i_RK,i_RK+1)*dt*Diffv];
%     Z  = [CD c_RK(i_RK)*dt*G; M Z1];
%     f  = [utemp; vtemp; -yM];
%     q  = Z\f;   
% 
%     uh = q(1:Nu);
%     vh = q(Nu+1:Nu+Nv);
%     p  = q(Nu+Nv+1:Nu+Nv+Np);
%     
%     
%     kp(:,i_RK) = p;
%     keyboard;

end

V      = [uh;vh];
p      = -3*kp(:,2) + 4*kp(:,3);
t      = tn;