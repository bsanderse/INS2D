% fully (energy) conservative time stepping
% not reversible in time (inviscid flow)

% like FC1, but splitting error is removed iteratively

CN     = 0.5;
f      = zeros(Nu+Nv+Np,1);
Z2     = spalloc(Np,Np,0);

% velocities at time level n
uhn        = uh;
vhn        = vh;
Vn         = V;
pn         = p;
tn         = t;

% estimate of the convecting velocity at n+1/2
% cu and cv should be the same for time n and time n+1 to conserve energy !!
if (EP==0)
    % Picard
    cu     = uh;
    cv     = vh;
elseif (EP==1)
    % extrapolated Picard:
    cu         = 1.5*uh - 0.5*uh_old;
    cv         = 1.5*vh - 0.5*vh_old;
end
% also possible for the pressure?

% evaluate BC and force at intermediate time 
t = tn + 0.5*dt;
boundary_conditions;
operator_boundary_conditions;
force;
t = tn;

error_nonlinear = 1;

% iterate to remove splitting error; fully conservative
% but not reversible, and maybe not 2nd order
% this could be done by solving the single linear(!) saddle-point system,
% but iterating until pressure correction error is removed is faster

i=1;
while (error_nonlinear > nonlinear_acc && i<nonlinear_maxit)
 

    
      utemp     = 0.5*(uh+uhn);
      vtemp     = 0.5*(vh+vhn);
      ptemp     = 0.5*(p+pn);

% uncommenting these lines will give reversibility in time (like time_FC2)
%       cu = utemp;
%       cv = vtemp;
    
      if (nonlinear_build_matrix == 1)
      % not completely sure of this implementation!
          error('not implemented'); 
      
      else
          CD_diag = spdiags(dt*Om_inv,0,Nu+Nv,Nu+Nv);
      end
      
      yConv_ux   = Cux*( (Iu_ux*cu+yIu_ux).*(Au_ux*utemp+yAu_ux) );
      yConv_uy   = Cuy*( (Iv_uy*cv+yIv_uy).*(Au_uy*utemp+yAu_uy) );
      yConv_vx   = Cvx*( (Iu_vx*cu+yIu_vx).*(Av_vx*vtemp+yAv_vx) );
      yConv_vy   = Cvy*( (Iv_vy*cv+yIv_vy).*(Av_vy*vtemp+yAv_vy) );
      
      % right-hand side, contains residual
      f(1:Nu)          = Omu.*(uh-uhn)/dt + yConv_ux + yConv_uy + ... 
                            -Diffu*utemp - yDiffu - Fx + Gx*ptemp + y_px;
      f(Nu+1:Nu+Nv)    = Omv.*(vh-vhn)/dt + yConv_vx + yConv_vy + ...
                            -Diffv*vtemp - yDiffv - Fy + Gy*ptemp + y_py;
      f(Nu+Nv+1:end)   = M*V + yM;
      

      % solve for dV and dp with Schur complement method
      % see Mullen et al. (2009)
      rhs = (-M*CD_diag*f(1:Nu+Nv) + f(Nu+Nv+1:end));
      
      if (nonlinear_build_matrix==0)
          % adapt rhs because in construction of A (see
          % operator_divergence) dt and CN are not present
          rhs = rhs/(dt*CN);
          % we can use precomputed CG or LU:
            if (poisson==1)
            % using pre-determined LU decomposition
              b   = L\rhs;                                                                                   
              dp  = U\b;
            elseif (poisson==2)
            % using preconditioned CG
              [dp, flag]  = pcg(-A,-rhs,CG_acc,CG_maxit,L,L',dp);
            elseif (poisson==3)  
              [dp,iter,norm1,norm2]=cg(B,int64(dia),int64(ndia),rhs,CG_acc,int64(Np),dp,int64(CG_maxit));
            elseif (poisson==4)
              [dp,iter,norm1,norm2]=cg_matlab(A,rhs,CG_acc,dp,CG_maxit,A_pc);
            end
      else
          A         = M*CD_diag*(CN*G);
          dp        = A\rhs;
      end
      
      dV        = CD_diag*(-f(1:Nu+Nv) - CN*G*dp);
      du        = dV(1:Nu);
      dv        = dV(Nu+1:Nu+Nv);

      
      V         = V + dV;
      uh        = uh + du;
      vh        = vh + dv;
      p         = p + dp;
      
      % divide the residual by a norm based on rhs
      error_nonlinear = norm(f)/norm(Om.*Vn/dt);
      
      i = i+1;
    
end

if (i>= nonlinear_maxit || error_nonlinear>nonlinear_acc || isnan(error_nonlinear))
    error(['nonlinear iterations did not converge in ' num2str(nonlinear_maxit) ' iterations']);
end

if (EP==1)
    uh_old    = uhn;
    vh_old    = vhn;
% p_old  = pn;
end
