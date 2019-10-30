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
%     
% we work with the following 'shifted' Butcher tableau, because A_RK(1,1)
% is always zero for explicit methods
A_RK = [A_RK(2:end,:); b_RK'];

% number of stages
s_RK = length(c_RK);

% store variables at start of time step
tn     = t;
uhn    = uh;
vhn    = vh;
pn     = p;

ku     = zeros(Nu,s_RK);
kv     = zeros(Nv,s_RK);

for i=1:s_RK
    
    t       = tn + c_RK(i)*dt;

    boundary_conditions;
    operator_boundary_conditions;
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
    utemp   = Omu_inv.*( - du2dx - duvdy + d2u + Fx);
    vtemp   = Omv_inv.*( - duvdx - dv2dy + d2v + Fy);
    R       = dt*[utemp;vtemp];

    % solve for the pressure (velocity update not used)
    pressure_correction;

    % velocity and pressure belonging to this stage
    p       = dp;

    % update conv+diffusion with pressure gradient
    ku(:,i) = utemp - Omu_inv.*(Gx*p + y_px);
    kv(:,i) = vtemp - Omv_inv.*(Gy*p + y_py);
    
    % update velocity current stage
    uh      = uhn + dt*ku*A_RK(i,:)';
    vh      = vhn + dt*kv*A_RK(i,:)';

end

t = tn; % time update is performed in solver_unsteady;

%% make pressure same order as velocity
if (p_add_solve==1)

    if (c_RK(end)~=1)
        t = tn + dtn;
        boundary_conditions;
        operator_boundary_conditions;
        force;
        t = tn;
    end
    
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