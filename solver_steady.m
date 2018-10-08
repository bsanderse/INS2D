% solve the entire saddlepoint system arising from the steady Navier-Stokes
% equations with linearization of the convective terms

if (restart.load == 0 && options.output.save_results)
    fprintf(fconv,'n            res              maxdiv           umom             vmom             k\n');
    fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e\n',...
        n,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
end

% zero block for total matrix
Np = options.grid.Np;
Nu = options.grid.Nu;
Nv = options.grid.Nv;

Z2 = spalloc(Np,Np,0);

% right hand side
f = zeros(Nu+Nv+Np,1);

Newton = 0;

maxres(2) = maxres(1);

G  = options.discretization.G;
M  = options.discretization.M;
yM = options.discretization.yM;

while ( maxres(n) > accuracy)
    
    % switch to Newton after nPicard steps
    if (strcmp(linearization,'Newton') && n>nPicard)
        Newton = 1;
    end
    
    n = n+1;
    
    %% formulate the solution of the nonlinear system in terms of the update
    % vector dq: 
    % dF dq = -res,
    % where dF is the Jacobian dFdq, and q contains both velocity and
    % pressure
    % res is the residual, which contains all the terms in the momentum
    % equation
    [maxres(n), fmom, dfmom] = F(uh,vh,p,t,options);
    fmass     = M*V+yM;
    f         = [fmom; fmass];
    Z         = [dfmom -G; M Z2];
    
    
    %% solve with direct solver from Matlab
    dq        = -Z\f;
    
    dV        = dq(1:Nu+Nv);
    dp        = dq(Nu+Nv+1:end);
    du        = dV(1:Nu);
    dv        = dV(Nu+1:end);
    
    V         = V + dV;
    uh        = uh + du;
    vh        = vh + dv;
    p         = p + dp;
    
    %% check residuals, conservation, write output files
    process_iteration;
    maxres(n)
    
    % write convergence information to file
    if (options.output.save_results == 1)
        fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e \n',...
            n,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
    end
end