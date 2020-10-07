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

while ( maxres(n) > options.solversettings.nonlinear_acc)
    
    % switch to Newton after nPicard steps if Jacobian_type=1
    if (options.solversettings.Jacobian_type == 1 && n>options.solversettings.nPicard)
        options.solversettings.Newton_factor = 1;
    end
    
    n = n+1;
    
    %% formulate the solution of the nonlinear system in terms of the update
    % vector dq: 
    % dF dq = -res,
    % where dF is the Jacobian dFdq, and q contains both velocity and
    % pressure
    % res is the residual, which contains all the terms in the momentum
    % equation and in the mass equations
    % fmom contains the right hand side of the momentum equation when
    % written in du/dt = fmom form, so fmom = -conv + diff - grad p
    [~, fmom, dfmom] = F(V,V,p,t,options,1);
    fmass     = M*V+yM;
    f         = [-fmom; fmass];
    % it is unclear why, but when using an asymmetric version of the G and
    % M blocks (using -M for mass equation), the Matlab solver converges much faster
    % 
    Z         = [dfmom -G; -M Z2];
    
    
    %% solve with direct solver
    dq        = Z\f;
    
    dV        = dq(1:Nu+Nv);
    dp        = dq(Nu+Nv+1:end);

    V         = V + dV;
    p         = p + dp;
    
    
    %% check residuals, conservation, write output files
    process_iteration;
    
    fprintf(fcw,'residual momentum equation: %16.8e \n', maxres(n));
    
    if (n>options.solversettings.nonlinear_maxit)
        warning(['Newton not converged in ' num2str(nonlinear_maxit) ' iterations, showing results anyway']);
        break
    end
    
end