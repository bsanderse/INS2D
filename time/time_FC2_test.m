% fully (energy) conservative time stepping
% solution of saddlepoint system,
% nonlinear_Newton iteration

% same as FC2, but more towards RK notation

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

utemp = uh;
vtemp = vh;

error_nonlinear = 1;
% iterate to remove linearization error; necessary to obtain full
% conservation and reversibility; if not removed still fully conservative
% but not reversible, and maybe not 2nd order
i=1;
while (error_nonlinear > nonlinear_acc  && i<nonlinear_maxit)
 
      % iterate until u^(n+1) = u^i = u^(i-1); this is different from
      % u^(n)!
    

%       utemp      = 0.5*(uh+uhn);
%       vtemp      = 0.5*(vh+vhn);
      ptemp      = 0.5*(p+pn);    
    
      %% convective terms, u-component
      % c^n * u^(n+1), c=u
      uIux       = Iu_ux*utemp+yIu_ux;                     % convective velocity, u_bar
      uAux       = Au_ux*utemp+yAu_ux;
      C1         = Cux*spdiags(uIux,0,N1,N1);   
      C2         = Cux*spdiags(uAux,0,N1,N1)*nonlinear_Newton;
      Conv_ux_11 = C1*Au_ux + C2*Iu_ux;
      yConv_ux   = Cux*( uIux.*uAux );

      vIuy       = Iv_uy*vtemp+yIv_uy;                     % convective velocity, v_bar
      uAuy       = Au_uy*utemp+yAu_uy;
      C1         = Cuy*spdiags(vIuy,0,N2,N2);   
      C2         = Cuy*spdiags(uAuy,0,N2,N2)*nonlinear_Newton;
      Conv_uy_11 = C1*Au_uy;
      Conv_uy_12 = C2*Iv_uy;
      yConv_uy   = Cuy*( vIuy.*uAuy );
            
      
      %% convective terms, v-component
      uIvx       = Iu_vx*utemp+yIu_vx;                 % convective velocity, u_bar  
      vAvx       = Av_vx*vtemp+yAv_vx;
      C1         = Cvx*spdiags(uIvx,0,N3,N3);
      C2         = Cvx*spdiags(vAvx,0,N3,N3)*nonlinear_Newton;
      Conv_vx_21 = C2*Iu_vx;
      Conv_vx_22 = C1*Av_vx;
      yConv_vx   = Cvx*( uIvx.*vAvx );
      
      vIvy       = Iv_vy*vtemp+yIv_vy;                 % convective velocity, v_bar
      vAvy       = Av_vy*vtemp+yAv_vy;
      C1         = Cvy*spdiags(vIvy,0,N4,N4);
      C2         = Cvy*spdiags(vAvy,0,N4,N4)*nonlinear_Newton;       
      Conv_vy_22 = C1*Av_vy + C2*Iv_vy;
      yConv_vy   = Cvy*( vIvy.*vAvy );
                         

      %% construct matrix (saddlepoint structure)
      CD_11      = spdiags(Omu,0,Nu,Nu) + ...
                   CN*dt*(- Diffu + Conv_ux_11 + Conv_uy_11);
      CD_12      = CN*dt*Conv_uy_12;
      CD_21      = CN*dt*Conv_vx_21;
      CD_22      = spdiags(Omv,0,Nv,Nv) + ...
                   CN*dt*(- Diffv + Conv_vx_22 + Conv_vy_22);    
      CD         = [CD_11 CD_12; CD_21 CD_22];
      Z          = [CD CN*G; M Z2]; 

      % force at midpoint
%       t  = t+0.5*dt;
%       force;
%       t  = tn;
      
      % right-hand side, contains residual
      f(1:Nu)          = Omu.*uh + yConv_ux + yConv_uy + ... 
                         -Diffu*utemp - yDiffu - Fx + Gx*ptemp + y_px;
      f(Nu+1:Nu+Nv)    = Omv.*vh + yConv_vx + yConv_vy  + ...
                         -Diffv*vtemp - yDiffv - Fy + Gy*ptemp + y_py;  
      f(Nu+Nv+1:end)   = M*V + yM;

      % solve
      dq        = -Z\f;
 
      dV        = dq(1:Nu+Nv);
      dp        = dq(Nu+Nv+1:end);
      du        = dV(1:Nu);
      dv        = dV(Nu+1:end);
      
      V         = V + dV;
      uh        = uh + du;
      vh        = vh + dv;
      p         = p + dp;

      utemp     = uhn + 0.5*dt*uh;
      vtemp     = vhn + 0.5*dt*vh;
      
      % divide the residual by a norm based on rhs
      error_nonlinear = norm(f)/norm(Om.*Vn/dt);
      
      i = i+1;
    
end
V = Vn + dt*V;
uh = V(1:Nu);
vh = V(Nu+1:end);

if (i>= nonlinear_maxit || error_nonlinear>nonlinear_acc || isnan(error_nonlinear))
    error(['nonlinear iterations did not converge in ' num2str(nonlinear_maxit) ' iterations']);
end