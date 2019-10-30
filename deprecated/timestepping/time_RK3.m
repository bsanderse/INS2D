% Runge-Kutta 3 without CN

% coefficients
beta1  = 4/15; beta2  =  1/15; beta3  = 1/6;
gamma1 = 0;    gamma2 =-17/60; gamma3 = -5/12;
sigma1 = 8/15; sigma2 =  5/12; sigma3 = 3/4;


dtn    = dt;

%% step 1
u1     = uh;
v1     = vh;
p1     = p;

% convection
cu     = u1;
cv     = v1;
convection;
conv_x1= du2dx + duvdy;
conv_y1= duvdx + dv2dy;

% diffusion
d2u1   = Diffu*u1 + yDiffu;
d2v1   = Diffv*v1 + yDiffv;


force;

u2     = u1 + dt*Omu_inv.*(- sigma1*(conv_x1 - d2u1) ...
                           - 2*beta1*(Gx*p1 + y_px) + sigma1*Fx);
v2     = v1 + dt*Omv_inv.*(- sigma1*(conv_y1 - d2v1) ...
                           - 2*beta1*(Gy*p1 + y_py) + sigma1*Fy);
                        
% u2     = (speye(Nu,Nu) - beta1*dt*spdiags(Omu_inv,0,Nu,Nu)*Diffu)\u2;
% v2     = (speye(Nv,Nv) - beta1*dt*spdiags(Omv_inv,0,Nv,Nv)*Diffv)\v2;

   
% enforce incompressibility 
R  = [u2;v2];
dt = 2*beta1*dtn;
pressure_correction;
u2 = V(1:Nu);
v2 = V(Nu+1:Nu+Nv);

p2 = p1 + dp;
dt = dtn;


%% step 2

% convection
uh     = u2;
vh     = v2;
cu     = u2;
cv     = v2;
convection;
conv_x2= du2dx + duvdy;
conv_y2= duvdx + dv2dy;

% diffusion
d2u2   = Diffu*u2 + yDiffu;
d2v2   = Diffv*v2 + yDiffv;

force;

u3     = u2 + dt*Omu_inv.*(- sigma2*(conv_x2 - d2u2) - gamma2*(conv_x1 - d2u1) ...
                           - 2*beta2*(Gx*p2 + y_px) + sigma2*Fx);
v3     = v2 + dt*Omv_inv.*(- sigma2*(conv_y2 - d2v2) - gamma2*(conv_y1 - d2v1) ...
                           - 2*beta2*(Gy*p2 + y_py) + sigma2*Fy);
                        
% u3     = (speye(Nu,Nu) - beta2*dt*spdiags(Omu_inv,0,Nu,Nu)*Diffu)\u3;
% v3     = (speye(Nv,Nv) - beta2*dt*spdiags(Omv_inv,0,Nv,Nv)*Diffv)\v3;

% enforce incompressibility 
R  = [u3;v3];
dt = 2*beta2*dtn;
pressure_correction;
u3 = V(1:Nu);
v3 = V(Nu+1:Nu+Nv);

p3 = p2 + dp;
dt = dtn;


%% step 3
% convection
uh     = u3;
vh     = v3;
cu     = u3;
cv     = v3;
convection;
conv_x3= du2dx + duvdy;
conv_y3= duvdx + dv2dy;

% diffusion
d2u3   = Diffu*u3 + yDiffu;
d2v3   = Diffv*v3 + yDiffv;

force;

uh     = u3 + dt*Omu_inv.*(- sigma3*(conv_x3 - d2u3) - gamma3*(conv_x2 - d2u2) ...
                            - 2*beta3*(Gx*p3 + y_px) + sigma3*Fx);
vh     = v3 + dt*Omv_inv.*(- sigma3*(conv_y3 - d2v3) - gamma3*(conv_y2 - d2v2) ...
                            - 2*beta3*(Gy*p3 + y_py) + sigma3*Fy);

% uh     = (speye(Nu,Nu) - beta3*dt*spdiags(Omu_inv,0,Nu,Nu)*Diffu)\uh;
% vh     = (speye(Nv,Nv) - beta3*dt*spdiags(Omv_inv,0,Nv,Nv)*Diffv)\vh;

% enforce incompressibility 
R  = [uh;vh];
dt = 2*beta3*dtn;
pressure_correction;
uh = V(1:Nu);
vh = V(Nu+1:Nu+Nv);
p  = p3 + dp;
dt = dtn;


if (p_add_solve==1)
    % additional poisson solve for 4th order accurate pressure
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