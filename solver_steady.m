% solve the entire saddlepoint system arising from the steady Navier-Stokes
% equations with linearization of the convective terms

if (restart.load == 0)
    fprintf(fconv,'n            res              maxdiv           umom             vmom             k\n');
    fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e\n',... 
                    n,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
end

% zero block for total matrix
Z2 = spalloc(Np,Np,0);

% right hand side
f = zeros(Nu+Nv+Np,1);

Newton = 0;

maxres(2) = maxres(1);

while ( maxres(n) > accuracy)

      % switch to Newton after nPicard steps
      if (strcmp(linearization,'Newton') && n>nPicard)
          Newton = 1;
      end

      n = n+1; 
      
      % to check if the convection matrix C is skew-symmetric, look at
%       Cu = Cux*spdiags(Iu_ux*uh+yIu_ux,0,N1,N1)*Au_ux + Cuy*spdiags(Iv_uy*vh+yIv_uy,0,N2,N2)*Au_uy;
%       Cv = Cvx*spdiags(Iu_vx*uh+yIu_vx,0,N3,N3)*Av_vx + Cvy*spdiags(Iv_vy*vh+yIv_vy,0,N4,N4)*Av_vy;
%       max2d(abs(Cu+Cu'))
%       max2d(abs(Cv+Cv'))


      %% convective terms, u-component
      % c^n * u^(n+1), c=u
      uIux       = Iu_ux*uh+yIu_ux;                     % convective velocity, u_bar
      uAux       = Au_ux*uh+yAu_ux;
      C1         = Cux*spdiags(uIux,0,N1,N1);   
      C2         = Cux*spdiags(uAux,0,N1,N1)*Newton;
      Conv_ux_11 = C1*Au_ux + C2*Iu_ux;
      yConv_ux   = C1*uAux;

      
      % c^n * u^(n+1), c=v
      vIuy       = Iv_uy*vh+yIv_uy;                     % convective velocity, v_bar
      uAuy       = Au_uy*uh+yAu_uy;
      C1         = Cuy*spdiags(vIuy,0,N2,N2);   
      C2         = Cuy*spdiags(uAuy,0,N2,N2)*Newton;
      Conv_uy_11 = C1*Au_uy;
      Conv_uy_12 = C2*Iv_uy;
      yConv_uy   = C1*uAuy;
      

      %% convective terms, v-component
      % c^n * v^(n+1), c=u
      uIvx       = Iu_vx*uh+yIu_vx;                 % convective velocity, u_bar  
      vAvx       = Av_vx*vh + yAv_vx;
      C1         = Cvx*spdiags(uIvx,0,N3,N3);
      C2         = Cvx*spdiags(vAvx,0,N3,N3)*Newton;
      Conv_vx_21 = C2*Iu_vx;
      Conv_vx_22 = C1*Av_vx;
      yConv_vx   = C1*vAvx;
      
      % c^n * v^(n+1), c=v
      vIvy       = Iv_vy*vh+yIv_vy;                 % convective velocity, v_bar
      vAvy       = Av_vy*vh+yAv_vy;
      C1         = Cvy*spdiags(vIvy,0,N4,N4);
      C2         = Cvy*spdiags(vAvy,0,N4,N4)*Newton;       
      Conv_vy_22 = C1*Av_vy + C2*Iv_vy;
      yConv_vy   = C1*vAvy;
      
      
      %% construct matrix (saddlepoint structure)
      CD_11      = - Diffu + Conv_ux_11 + Conv_uy_11;
      CD_12      = Conv_uy_12;
      CD_21      = Conv_vx_21;
      CD_22      = - Diffv + Conv_vx_22 + Conv_vy_22;     
      CD         = [CD_11 CD_12; CD_21 CD_22];
      Z          = [CD G; M Z2];
      
      Z          = Z + relax*speye(Nu+Nv+Np,Nu+Nv+Np);
      
      % right-hand side; this is -1*residual
      f(1:Nu)        = Diffu*uh + yDiffu - yConv_ux - yConv_uy + ...
                       Fx - Gx*p - y_px;
      f(Nu+1:Nu+Nv)  = Diffv*vh + yDiffv - yConv_vx - yConv_vy + ...
                       Fy - Gy*p - y_py;
      f(Nu+Nv+1:end) = -M*V - yM;
      
      
      % using Schur complement (too slow for steady problems)
%       CD_diag   = spdiags(1./diag(CD),0,Nu+Nv,Nu+Nv);   
%       Z         = M*CD_diag*G;
%       dp        = Z\(M*CD_diag*f(1:Nu+Nv) - f(Nu+Nv+1:end));
%       dV        = CD_diag*(f(1:Nu+Nv) - G*dp);
%       du        = dV(1:Nu);
%       dv        = dV(Nu+1:end);
      
      %% solve with direct solver from Matlab
      dq        = Z\f;
     
      dV        = dq(1:Nu+Nv);
      dp        = dq(Nu+Nv+1:end);
      du        = dV(1:Nu);
      dv        = dV(Nu+1:end);     
            
      V         = V + dV;
      uh        = uh + du;
      vh        = vh + dv;
      p         = p + dp;

      % check residuals, conservation, write output files    
      process_iteration;
      maxres(n)

      % write convergence information to file  
      fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e \n',...
                    n,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
%   pause;
end