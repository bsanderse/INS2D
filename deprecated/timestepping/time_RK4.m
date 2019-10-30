% Runge-Kutta 4
b_RK4 = [1;2;2;1]/6;
% B_RK4 = diag(b_RK4);
a_RK4 = [1; 1; 2]/2;
A_RK4 = spdiags(a_RK4,-1,4,4);
% M_RK4 = (B_RK4*A_RK4+A_RK4'*B_RK4-b_RK4*b_RK4');
c_RK4 = [0; 1/2; 1/2; 1];

% store variables at start of time step
tn     = t;
uhn    = uh;
vhn    = vh;
pn     = p;

dtn    = dt;


%% step 1
t      = tn;
uh     = uhn;
vh     = vhn;

boundary_conditions;
operator_boundary_conditions;

% convection
cu     = uh;
cv     = vh;
convection;

% diffusion
d2u    = Diffu*uh + yDiffu;
d2v    = Diffv*vh + yDiffv;

force;

ku1    = Omu_inv.*( - du2dx - duvdy + d2u + Fx - Gx*p - y_px);
kv1    = Omv_inv.*( - duvdx - dv2dy + d2v + Fy - Gy*p - y_py);
p1     = p;


%% step 2
dt2    = c_RK4(2)*dtn;
t      = tn  + dt2;
uh     = uhn + dt*A_RK4(2,1)*ku1;    
vh     = vhn + dt*A_RK4(2,1)*kv1;

boundary_conditions;
operator_boundary_conditions;


if (method==82)
    R  = [uh;vh];
    % enforce incompressibility 
    dt = dt2;
    pressure_correction;
    uh = V(1:Nu);
    vh = V(Nu+1:Nu+Nv);
    % Kampanis:
%     p  = p1 + dp;    
    p  = pn + dp;
    dt = dtn;
    % returns uh, vh, p
end

uh2    = uh;
vh2    = vh;

% convection
cu     = uh;
cv     = vh;
convection;

% diffusion
d2u    = Diffu*uh + yDiffu;
d2v    = Diffv*vh + yDiffv;

force;

ku2    = Omu_inv.*( - du2dx - duvdy + d2u + Fx - Gx*p - y_px);
kv2    = Omv_inv.*( - duvdx - dv2dy + d2v + Fy - Gy*p - y_py);
p2     = p;


%% step 3
dt3    = c_RK4(3)*dtn;
t      = tn  + dt3;
uh     = uhn + dt*A_RK4(3,2)*ku2;    
vh     = vhn + dt*A_RK4(3,2)*kv2;

boundary_conditions;
operator_boundary_conditions;

if (method==82)
    R  = [uh;vh];
    % enforce incompressibility 
    dt = dt3;    
    pressure_correction;
    uh = V(1:Nu);
    vh = V(Nu+1:Nu+Nv);
    % Kampanis:
%     p  = p2 + dp;
    p  = pn + dp;
    dt = dtn;
    % returns uh, vh, p
end

uh3    = uh;
vh3    = vh;


% convection
cu     = uh;
cv     = vh;
convection;

% diffusion
d2u    = Diffu*uh + yDiffu;
d2v    = Diffv*vh + yDiffv;

ku3    = Omu_inv.*( - du2dx - duvdy + d2u + Fx - Gx*p - y_px);
kv3    = Omv_inv.*( - duvdx - dv2dy + d2v + Fy - Gy*p - y_py);
p3     = p;


%% step 4
dt4    = c_RK4(4)*dt;
t      = tn  + dt4;
uh     = uhn + dt*A_RK4(4,3)*ku3;    
vh     = vhn + dt*A_RK4(4,3)*kv3;

boundary_conditions;
operator_boundary_conditions;

if (method==82)
    R  = [uh;vh];
    % enforce incompressibility 
    dt = dt4;    
    pressure_correction;
    uh = V(1:Nu);
    vh = V(Nu+1:Nu+Nv);
    % Kampanis:
%     p  = p3 + dp;
    p  = pn + dp;
    dt = dtn;
    % returns uh, vh, p
end

uh4    = uh;
vh4    = vh;

% convection
cu     = uh;
cv     = vh;
convection;

% diffusion
d2u    = Diffu*uh + yDiffu;
d2v    = Diffv*vh + yDiffv;

ku4    = Omu_inv.*( - du2dx - duvdy + d2u + Fx - Gx*p - y_px);
kv4    = Omv_inv.*( - duvdx - dv2dy + d2v + Fy - Gy*p - y_py);
p4     = p;


%% final step
uh     = uhn + dt*(b_RK4(1)*ku1 + b_RK4(2)*ku2 + b_RK4(3)*ku3 + b_RK4(4)*ku4);
vh     = vhn + dt*(b_RK4(1)*kv1 + b_RK4(2)*kv2 + b_RK4(3)*kv3 + b_RK4(4)*kv4);

R      = [uh;vh];

t      = tn; % back to old time step, will be updated in solver_unsteady
dt     = dtn;

pressure_correction;

% Kampanis:
% pK     = (1/6)*(p1 + 2*p2 + 2*p3 + p4);
% p      = pK + dp;

p      = pn + dp;

uh     = V(1:Nu);
vh     = V(Nu+1:Nu+Nv);


%% evaluate energy error according to theory
% 
% % Ainv   = A\speye(Np);
% Ainv   = (U\speye(Np))*(L\speye(Np));
% Q      = speye(Nu+Nv) - spdiags(Om_inv,0,Nu+Nv,Nu+Nv)*G*Ainv*M;
% k1     = [ku1;kv1]; k2 = [ku2;kv2]; k3 = [ku3;kv3]; k4 = [ku4;kv4];
% V1     = [uh1;vh1]; V2 = [uh2;vh2]; V3 = [uh3;vh3]; V4 = [uh4;vh4];
% Vhn    = Q*[uhn;vhn]; % this has no or little effect, because Vhn should be div-free
% k1     = Q*k1; k2 = Q*k2; k3 = Q*k3; k4 = Q*k4;
% 
% k(n)   = 0.5*sum(Omu.*uh.^2) + 0.5*sum(Omv.*vh.^2);
% k(n) - k(n-1)
% 
% 
% k_2 = M_RK4(1,1)*k1.^2 + M_RK4(2,2)*k2.^2 + ...
%       M_RK4(3,3)*k3.^2 + M_RK4(4,4)*k4.^2 + ...
%       M_RK4(1,2)*2*k1.*k2 + M_RK4(1,3)*2*k1.*k3 + ...
%       M_RK4(1,4)*2*k1.*k4 + M_RK4(2,3)*2*k2.*k3 + ...
%       M_RK4(2,4)*2*k2.*k4 + M_RK4(3,4)*2*k3.*k4;                      
%                       
% 
% 0.5*sum(Om.*(Vhn).^2) + 0.5*2*dt*sum( Om.*(b_RK4(1)*k1.*V1 + b_RK4(2)*k2.*V2 + ...
%                                            b_RK4(3)*k3.*V3 + b_RK4(4)*k4.*V4) ) + ...
%                       - 0.5*(dt^2)*sum(Om.*k_2) - k(n-1)
% 
% 
% keyboard;
% 

% we should have sum(kui.*uhi)+sum(kvi.*vhi)=0 for each i


%%
if (p_add_solve==1)
    % force and boundary conditions need not to be updated, because
    % last stage has c_4=1, so the BC are at the right time
%     t = tn + dtn;
%     force;
%     boundary_conditions;
%     operator_boundary_conditions;
%     t = tn;
    
    % additional poisson solve for 4th order accurate pressure
    % convection
    cu = uh;
    cv = vh;
    convection;

    % diffusion
    d2u    = Diffu*uh + yDiffu;
    d2v    = Diffv*vh + yDiffv;


    Ru =  - du2dx - duvdy + d2u + Fx - y_px;    
    Rv =  - duvdx - dv2dy + d2v + Fy - y_py;
    R  = Om_inv.*[Ru;Rv];
    f  = (M*R + yM);

    dp = zeros(Np,1);

    if (poisson==1)
        % using pre-determined LU decomposition
        b   = L\f;                                                                                   
        dp  = U\b;
    elseif (poisson==2)
      % using preconditioned CG
        t_pressure = toc;
        [dp, flag, norm1, iter]  = pcg(-A,-f,CG_acc,CG_maxit,L,L',dp);
        t_solve    = toc - t_pressure;
    elseif (poisson==3)  
        t_pressure = toc;
        [dp,iter,norm1,norm2]=cg(B,int64(dia),int64(ndia),f,CG_acc,int64(Np),dp,int64(CG_maxit));
        t_solve    = toc - t_pressure;    
    elseif (poisson==4)
        t_pressure = toc;  
        [dp,iter,norm1,norm2]=cg_matlab(A,f,CG_acc,dp,CG_maxit,A_pc);
        t_solve    = toc - t_pressure;
    elseif (poisson==5)
        t_pressure = toc;
        PetscBinaryWrite(PS,-f');
        t_write    = toc - t_pressure;

        t_pressure = toc;   
        dp     = PetscBinaryRead(PS);
        iter   = PetscBinaryRead(PS);
        norm1  = PetscBinaryRead(PS);   
        t_read = toc - t_pressure;

        dp     = dp';
        iter   = iter(1);
        norm1  = norm1(1);  
    end

    % write information on pressure solve
    if (poisson==2 || poisson==3 || poisson==4)
        fprintf(fpres,'%-10i %-10i %16.8e %16.8e\n',n,iter,norm1,t_solve);
    elseif (poisson == 5)
        fprintf(fpres,'%-10i %-10i %16.8e %16.8e %16.8e\n',n,iter,norm1,t_write,t_read);
    end    

    p = dp;
end