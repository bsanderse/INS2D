% fully (energy) conservative time stepping
% solution of saddlepoint system,
% nonlinear_Newton iteration

% reversible in time (inviscid flow)

CN     = 0.5;
f      = zeros(Nu+Nv+Np,1);
Z2     = spalloc(Np,Np,0);

% velocities at time level n
uhn    = uh;
vhn    = vh;
Vn     = V;
pn     = p;
tn     = t;

error_nonlinear = 1;
% iterate to remove linearization error; necessary to obtain full
% conservation and reversibility; if not removed still fully conservative
% but not reversible, and maybe not 2nd order
i=1;

% evaluate BC and force for momentum equation at intermediate time 
t = tn + 0.5*dt;

boundary_conditions;
interpolate_bc;
operator_bc_momentum;
force;

% evaluate BC for divergence at next time level
t = tn + dt;
boundary_conditions;
interpolate_bc;
operator_bc_divergence;

t = tn;

while (error_nonlinear > nonlinear_acc  && i<nonlinear_maxit)
 
      % iterate until u^(n+1) = u^i = u^(i-1); this is different from
      % u^(n)!
    
      utemp      = 0.5*(uh+uhn);
      vtemp      = 0.5*(vh+vhn);
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
      CD_11      = spdiags(Omu,0,Nu,Nu)/dt + ...
                   CN*(- Diffu + Conv_ux_11 + Conv_uy_11);
      CD_12      = CN*Conv_uy_12;
      CD_21      = CN*Conv_vx_21;
      CD_22      = spdiags(Omv,0,Nv,Nv)/dt + ...
                   CN*(- Diffv + Conv_vx_22 + Conv_vy_22);    
      CD         = [CD_11 CD_12; CD_21 CD_22];
      Z          = [CD CN*G; M Z2]; 
      
      % right-hand side, contains residual
      f(1:Nu)          = Omu.*(uh-uhn)/dt + yConv_ux + yConv_uy + ... 
                         -Diffu*utemp - yDiffu - Fx + Gx*ptemp + y_px;
      f(Nu+1:Nu+Nv)    = Omv.*(vh-vhn)/dt + yConv_vx + yConv_vy  + ...
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

      % note that dp does not necessarily go to zero
      
      % divide the residual by a norm based on rhs
%       error_nonlinear = norm(f)/norm(Om.*Vn/dt);
      error_nonlinear = norm(f,inf);
%       error_nonlinear = norm(dq)
      i = i+1;
    
end

if (i>= nonlinear_maxit || error_nonlinear>nonlinear_acc || isnan(error_nonlinear))
    error(['nonlinear iterations did not converge in ' num2str(nonlinear_maxit) ' iterations']);
end

if (p_add_solve==1)
    % evaluate BC and force at final time    
    t = tn + dt;
    boundary_conditions;
    interpolate_bc;
    operator_bc_momentum;
    operator_bc_divergence;
    force;
    t = tn;    
    
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
    f  = M*R + ydM;

    pressure_poisson;

    p = dp;
end
