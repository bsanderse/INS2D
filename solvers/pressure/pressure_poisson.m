function dp = pressure_poisson(f,t,options)
% compute pressure from pressure poisson problem with right-hand side f

% assume the Laplace matrix is known (A) and is possibly factorized (LU);
% right hand side is given by f
% we should have sum(f) = 0 for periodic and no-slip BC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Solve the Poisson equation for the pressure                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np   = options.grid.Np;
dp   = zeros(Np,1);

fpres = options.output.fpres;

poisson  = options.solversettings.poisson;
CG_acc   = options.solversettings.CG_acc;
CG_maxit = options.solversettings.CG_maxit;

% check if a Poisson solve is necessary
A   = options.discretization.A;
tol = 1e-14;
if (max(abs(A*dp-f))<tol)
    return;
end

% choose a method to solve the Poisson equation
if (poisson==1)
    % using pre-determined LU decomposition
    L   = options.solversettings.L;
    U   = options.solversettings.U;
    b   = L\f;
    dp  = U\b;
elseif (poisson==2)
    % using preconditioned CG
    L   = options.solversettings.L;
    
    t_pressure = toc;
    [dp, flag, norm1, iter]  = pcg(-A,-f,CG_acc,CG_maxit,L,L',dp);
    if (flag>0)
        warning('PCG not converged');
        if (flag~=3)
            keyboard
        end
    end
    t_solve    = toc - t_pressure;
elseif (poisson==3)
    B    = options.solversettings.B;
    dia  = options.solversettings.dia;
    ndia = options.solversettings.ndia;
    t_pressure = toc;
    [dp,iter,norm1,norm2]=cg(B,int64(dia),int64(ndia),f,CG_acc,int64(Np),dp,int64(CG_maxit));
    t_solve    = toc - t_pressure;
elseif (poisson==4)
    A_pc   = options.solversettings.A_pc;
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
    %    t_solve= PetscBinaryRead(PS);
    t_read = toc - t_pressure;
    iter   = iter(1);
    norm1  = norm1(1);
    
    dp     = dp';
    
else
    error('wrong poisson method');
end

% write information on pressure solve
if (options.output.save_results == 1)
    if (poisson==2 || poisson==3 || poisson==4)
        fprintf(fpres,'%-10i %-10i %16.8e %16.8e\n',t,iter,norm1,t_solve);
    elseif (poisson == 5)
        fprintf(fpres,'%-10i %-10i %16.8e %16.8e %16.8e\n',t,iter,norm1,t_write,t_read);
    end
end

end