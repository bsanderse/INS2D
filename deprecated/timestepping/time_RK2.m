% Runge-Kutta 2, Heun
b_RK2 = [1;1]/2;
B_RK2 = diag(b_RK2);
a_RK2 = [1];
A_RK2 = spdiags(a_RK2,-1,2,2);
M_RK2 = -(B_RK2*A_RK2+A_RK2'*B_RK2-b_RK2*b_RK2');


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
t      = tn  + dt;
uh     = uhn + dt*ku1;    
vh     = vhn + dt*kv1;

boundary_conditions;
operator_boundary_conditions;

if (method==102)
    R  = [uh;vh];
    % enforce incompressibility 
    dt = dtn;
    pressure_correction;
    uh = V(1:Nu);
    vh = V(Nu+1:Nu+Nv);
    p  = pn + dp;
    dt = dtn;
    % returns uh, vh, p
end

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


%% final step
uh     = uhn + (1/2)*dt*(ku1 + ku2);
vh     = vhn + (1/2)*dt*(kv1 + kv2);

R      = [uh;vh];

t      = tn; % back to old time step, will be updated in solver_unsteady
dt     = dtn;
   
pressure_correction;

p      = pn + dp;

uh     = V(1:Nu);
vh     = V(Nu+1:Nu+Nv);


%% evaluate energy error according to theory
% we should have sum(kui.*uhi)+sum(kvi.*vhi)=0 for each i



% ku and kv are divided by volume
k_error_u = M_RK2(1,1)*ku1.*ku1 + M_RK2(2,2)*ku2.*ku2 + ...
            M_RK2(1,2)*ku1.*ku2 + M_RK2(2,1)*ku2.*ku1;
k_error_v = M_RK2(1,1)*kv1.*kv1 + M_RK2(2,2)*kv2.*kv2 + ...
            M_RK2(1,2)*kv1.*kv2 + M_RK2(2,1)*kv2.*kv1;
        
dk_error(n) = (sum(Omu.*k_error_u) + sum(Omv.*k_error_v))*(dt^2)*0.5;



if (p_add_solve==1)
    % force and boundary conditions need not to be updated, because
    % last stage has c_4=1, so the BC are at the right time
%     t = tn + dtn;
%     force;
%     boundary_conditions;
%     operator_boundary_conditions;
%     t = tn;
        
    
    % additional poisson solve for 2nd order accurate pressure
    % convection
    cu = uh;
    cv = vh;
    convection;

    % diffusion
    d2u    = Diffu*uh + yDiffu;
    d2v    = Diffv*vh + yDiffv;

    % force from current time level
    tn = t;
    t  = t+dt;
    force;
    t  = tn;

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
        [dp,iter,norm1,norm2]=cg_matlab_gen(A,f,CG_acc,Np,dp,CG_maxit,A_pc);
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