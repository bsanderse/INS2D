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
      Cu = Cux*spdiags(Iu_ux*uh+yIu_ux,0,N1,N1)*Au_ux + Cuy*spdiags(Iv_uy*vh+yIv_uy,0,N2,N2)*Au_uy;
      Cv = Cvx*spdiags(Iu_vx*uh+yIu_vx,0,N3,N3)*Av_vx + Cvy*spdiags(Iv_vy*vh+yIv_vy,0,N4,N4)*Av_vy;
      Cu3 = Cux3*spdiags(Iu_ux3*uh+yIu_ux3,0,length(yIu_ux3),length(yIu_ux3))*Au_ux3 + Cuy3*spdiags(Iv_uy3*vh+yIv_uy3,0,length(yIv_uy3),length(yIv_uy3))*Au_uy3;
      Cv3 = Cvx3*spdiags(Iu_vx3*uh+yIu_vx3,0,length(yIu_vx3),length(yIu_vx3))*Av_vx3 + Cvy3*spdiags(Iv_vy3*vh+yIv_vy3,0,length(yIv_vy3),length(yIv_vy3))*Av_vy3; 
% 
% %       then alfa*Cu - Cu3 and alfa*Cv - Cv3 are skew-symmetric     
%       max2d(abs( (alfa*Cu-Cu3)' + (alfa*Cu-Cu3)))
%       max2d(abs( (alfa*Cv-Cv3)' + (alfa*Cv-Cv3)))
% keyboard
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
      
      
      if (order4==1)
          %% u-component
          uIux3  = Iu_ux3*uh+yIu_ux3;
          uAux3  = Au_ux3*uh+yAu_ux3;
          C1     = Cux3*spdiags(uIux3,0,length(uIux3),length(uIux3));   
          C2     = Cux3*spdiags(uAux3,0,length(uAux3),length(uAux3))*Newton;
          Conv_ux_11_3 = C1*Au_ux3 + C2*Iu_ux3;
          yConv_ux_3   = C1*uAux3;
          
          vIuy3  = Iv_uy3*vh+yIv_uy3;
          uAuy3  = Au_uy3*uh+yAu_uy3;
          C1         = Cuy3*spdiags(vIuy3,0,length(vIuy3),length(vIuy3));   
          C2         = Cuy3*spdiags(uAuy3,0,length(vIuy3),length(vIuy3))*Newton;
          Conv_uy_11_3 = C1*Au_uy3;
          Conv_uy_12_3 = C2*Iv_uy3;
          yConv_uy_3   = C1*uAuy3;          
          
          %%  v-component
          uIvx3      = Iu_vx3*uh+yIu_vx3;                 % convective velocity, u_bar  
          vAvx3      = Av_vx3*vh+yAv_vx3;
          C1         = Cvx3*spdiags(uIvx3,0,length(uIvx3),length(uIvx3));
          C2         = Cvx3*spdiags(vAvx3,0,length(vAvx3),length(vAvx3))*Newton;
          Conv_vx_21_3 = C2*Iu_vx3;
          Conv_vx_22_3 = C1*Av_vx3;
          yConv_vx_3   = C1*vAvx3;

          vIvy3      = Iv_vy3*vh+yIv_vy3;                 % convective velocity, v_bar
          vAvy3      = Av_vy3*vh+yAv_vy3;
          C1         = Cvy3*spdiags(vIvy3,0,length(vIvy3),length(vIvy3));
          C2         = Cvy3*spdiags(vAvy3,0,length(vAvy3),length(vAvy3))*Newton;       
          Conv_vy_22_3 = C1*Av_vy3 + C2*Iv_vy3;
          yConv_vy_3   = C1*vAvy3;
      
      end
      
      %% construct matrix (saddlepoint structure)
      
      if (order4==0)
          CD_11  = - Diffu + Conv_ux_11 + Conv_uy_11;
          CD_12  = Conv_uy_12;
          CD_21  = Conv_vx_21;
          CD_22  = - Diffv + Conv_vx_22 + Conv_vy_22;     
      elseif (order4==1)
          CD_11  = - Diffu + alfa*Conv_ux_11 - Conv_ux_11_3 + alfa*Conv_uy_11 - Conv_uy_11_3;
          CD_12  = alfa*Conv_uy_12 - Conv_uy_12_3;
          CD_21  = alfa*Conv_vx_21 - Conv_vx_21_3;
          CD_22  = - Diffv + alfa*Conv_vx_22 - Conv_vx_22_3 + alfa*Conv_vy_22 - Conv_vy_22_3;
      end
      
      CD         = [CD_11 CD_12; CD_21 CD_22];
      Z          = [CD G; M Z2];
      
      Z          = Z + relax*speye(Nu+Nv+Np,Nu+Nv+Np);
      
      % right-hand side; this is -1*residual
      if (order4==0)
          f(1:Nu)        = Diffu*uh + yDiffu - yConv_ux - yConv_uy + ...
                           Fx - Gx*p - y_px;
          f(Nu+1:Nu+Nv)  = Diffv*vh + yDiffv - yConv_vx - yConv_vy + ...
                           Fy - Gy*p - y_py;
      elseif (order4==1)
          f(1:Nu)        = Diffu*uh + yDiffu + ...
                           - (alfa*yConv_ux - yConv_ux_3) + ...
                           - (alfa*yConv_uy - yConv_uy_3) + ...
                           Fx - Gx*p - y_px;
          f(Nu+1:Nu+Nv)  = Diffv*vh + yDiffv + ...
                           - (alfa*yConv_vx - yConv_vx_3) + ...
                           - (alfa*yConv_vy - yConv_vy_3) + ...
                           Fy - Gy*p - y_py;         
      end
      f(Nu+Nv+1:end) = -M*V - yM;

      
      % using Schur complement (too slow for steady problems)
%       CD_diag   = spdiags(1./diag(CD),0,Nu+Nv,Nu+Nv);   
%       Z         = M*CD_diag*G;
%       dp        = Z\(M*CD_diag*f(1:Nu+Nv) - f(Nu+Nv+1:end));
%       dV        = CD_diag*(f(1:Nu+Nv) - G*dp);
%       du        = dV(1:Nu);
%       dv        = dV(Nu+1:end);
% 
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
%       max(abs(f))
% keyboard
      % write convergence information to file  
      fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e \n',...
                    n,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
%   pause;
end