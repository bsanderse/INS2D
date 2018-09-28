% assume the Laplace matrix is known (A) and is possibly factorized (LU);
% right hand side is given by f
% we should have sum(f) = 0 for periodic and no-slip BC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Solve the Poisson equation for the pressure                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
dp   = zeros(Np,1);


if (poisson==1) 
    % using pre-determined LU decomposition
    b   = L\f;                                                                                   
    dp  = U\b;
elseif (poisson==2)
    % using preconditioned CG
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
%    t_solve= PetscBinaryRead(PS);
    t_read = toc - t_pressure;
    iter   = iter(1);
    norm1  = norm1(1);
    
    dp     = dp';

else
    error('wrong poisson method');
end

% write information on pressure solve
if (poisson==2 || poisson==3 || poisson==4)
    fprintf(fpres,'%-10i %-10i %16.8e %16.8e\n',n,iter,norm1,t_solve);
elseif (poisson == 5)
    fprintf(fpres,'%-10i %-10i %16.8e %16.8e %16.8e\n',n,iter,norm1,t_write,t_read);
end