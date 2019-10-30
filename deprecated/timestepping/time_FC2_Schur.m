% fully (energy) conservative time stepping
% solution of saddlepoint system,
% nonlinear_Newton iteration

% reversible in time (inviscid flow)

CN     = 0.5;
f      = zeros(Nu+Nv+Np,1);
Z2     = spalloc(Np,Np,0);

% velocities at time level n
uhn        = uh;
vhn        = vh;
Vn         = V;
pn         = p;
tn         = t;

% evaluate BC and force at intermediate time 
t = tn + 0.5*dt;
boundary_conditions;
operator_boundary_conditions;
force;
t = tn;

error_nonlinear = 1;

% iterate to remove linearization error; necessary to obtain full
% conservation and reversibility; if not removed still fully conservative
% but not reversible, and maybe not 2nd order
i=1;
while (error_nonlinear > nonlinear_acc  && i<nonlinear_maxit)
 
      % iterate until u^(n+1) = u^i = u^(i-1); this is different from
      % u^(n)!
    
      utemp     = 0.5*(uh+uhn);
      vtemp     = 0.5*(vh+vhn);
      ptemp     = 0.5*(p+pn);    
    
      if (nonlinear_build_matrix == 1)
      % not completely sure of this implementation!      
      
          %% convective terms, u-component
          % c^n * u^(n+1), c=u
          uIux       = Iu_ux*utemp+yIu_ux;                     % convective velocity, u_bar
          uAux       = Au_ux*utemp+yAu_ux;
          C1         = Cux*spdiags(uIux,0,N1,N1);   
          C2         = Cux*spdiags(uAux,0,N1,N1)*nonlinear_Newton;
          Conv_ux_11 = C1*Au_ux + C2*Iu_ux;

          vIuy       = Iv_uy*vtemp+yIv_uy;                     % convective velocity, v_bar
          uAuy       = Au_uy*utemp+yAu_uy;
          C1         = Cuy*spdiags(vIuy,0,N2,N2);   
          C2         = Cuy*spdiags(uAuy,0,N2,N2)*nonlinear_Newton;
          Conv_uy_11 = C1*Au_uy;
          Conv_uy_12 = C2*Iv_uy;


          %% convective terms, v-component
          uIvx       = Iu_vx*utemp+yIu_vx;                 % convective velocity, u_bar  
          vAvx       = Av_vx*vtemp+yAv_vx;
          C1         = Cvx*spdiags(uIvx,0,N3,N3);
          C2         = Cvx*spdiags(vAvx,0,N3,N3)*nonlinear_Newton;
          Conv_vx_21 = C2*Iu_vx;
          Conv_vx_22 = C1*Av_vx;

          vIvy       = Iv_vy*vtemp+yIv_vy;                 % convective velocity, v_bar
          vAvy       = Av_vy*vtemp+yAv_vy;
          C1         = Cvy*spdiags(vIvy,0,N4,N4);
          C2         = Cvy*spdiags(vAvy,0,N4,N4)*nonlinear_Newton;       
          Conv_vy_22 = C1*Av_vy + C2*Iv_vy;



          %% construct matrix (saddlepoint structure)
          CD_11      = spdiags(Omu,0,Nu,Nu)/dt + ...
                       CN*(- Diffu + Conv_ux_11 + Conv_uy_11);
          CD_12      = CN*Conv_uy_12;
          CD_21      = CN*Conv_vx_21;
          CD_22      = spdiags(Omv,0,Nv,Nv)/dt + ...
                       CN*(- Diffv + Conv_vx_22 + Conv_vy_22);    
          CD         = [CD_11 CD_12; CD_21 CD_22];

          CD_diag    = spdiags(1./diag(CD),0,Nu+Nv,Nu+Nv);   

      else
           CD_diag = spdiags(dt*Om_inv,0,Nu+Nv,Nu+Nv);
      end
      
      yConv_ux   = Cux*( (Iu_ux*utemp+yIu_ux).*(Au_ux*utemp+yAu_ux) );
      yConv_uy   = Cuy*( (Iv_uy*vtemp+yIv_uy).*(Au_uy*utemp+yAu_uy) );
      yConv_vx   = Cvx*( (Iu_vx*utemp+yIu_vx).*(Av_vx*vtemp+yAv_vx) );
      yConv_vy   = Cvy*( (Iv_vy*vtemp+yIv_vy).*(Av_vy*vtemp+yAv_vy) );

      % right-hand side, contains residual
      f(1:Nu)          = Omu.*(uh-uhn)/dt + yConv_ux + yConv_uy + ... 
                         -Diffu*utemp - yDiffu - Fx + Gx*ptemp + y_px;
      f(Nu+1:Nu+Nv)    = Omv.*(vh-vhn)/dt + yConv_vx + yConv_vy  + ...
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
              [dp, flag, relres, iter]  = pcg(-A,-rhs,CG_acc,CG_maxit,L,L',dp);
              if (flag>0)
                  error('PCG not converged');
              end
            %            [x, options] = AMGsolver(-A, PREC, options, -b, dp); 
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
      dv        = dV(Nu+1:end);
      
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