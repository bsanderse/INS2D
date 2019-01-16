% Implicit Midpoint with iteration to remove linearization error
% nonlinear_Newton iteration
% pressure correction

CN     = 0.5;
f      = zeros(Nu+Nv,1);

% velocities at old time level n
uhn    = uh;
vhn    = vh;
Vn     = V;
tn     = t;
pn     = p;

% ptemp      = 1.5*p - 0.5*p_old;
ptemp = p;

error_nonlinear = 1;
% iterate to remove linearization error (still splitting error from
% pressure correction present)
i=1;

% evaluate BC and force at midpoint
t = tn + 0.5*dt;

boundary_conditions;
interpolate_bc;
% operator_bc_divergence;
operator_bc_momentum;
force;

t = tn;

while (error_nonlinear > nonlinear_acc && i<nonlinear_maxit)
 
      % iterate until u^(n+1) = u^i = u^(i-1); this is different from
      % u^(n)!
    
      utemp     = 0.5*(uh+uhn);
      vtemp     = 0.5*(vh+vhn);
       
      if (nonlinear_build_matrix == 1)

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
                       CN*(-Diffu + Conv_ux_11 + Conv_uy_11); %- Diffu 
          CD_12      = CN*Conv_uy_12;
          CD_21      = CN*Conv_vx_21;
          CD_22      = spdiags(Omv,0,Nv,Nv)/dt + ...
                       CN*(-Diffv + Conv_vx_22 + Conv_vy_22); %- Diffv   
          CD         = [CD_11 CD_12; CD_21 CD_22];
      else
          CD         = spdiags(Om/dt,0,Nu+Nv,Nu+Nv);
      end
      
      yConv_ux   = Cux*( (Iu_ux*utemp+yIu_ux).*(Au_ux*utemp+yAu_ux) );
      yConv_uy   = Cuy*( (Iv_uy*vtemp+yIv_uy).*(Au_uy*utemp+yAu_uy) );
      yConv_vx   = Cvx*( (Iu_vx*utemp+yIu_vx).*(Av_vx*vtemp+yAv_vx) );
      yConv_vy   = Cvy*( (Iv_vy*vtemp+yIv_vy).*(Av_vy*vtemp+yAv_vy) );
            
      % right-hand side, contains residual
      f(1:Nu)             = Omu.*(uh-uhn)/dt + yConv_ux + yConv_uy +... 
                            -Diffu*utemp - yDiffu - Fx + Gx*ptemp + y_px;
      f(Nu+1:Nu+Nv)       = Omv.*(vh-vhn)/dt + yConv_vx + yConv_vy + ...
                            -Diffv*vtemp - yDiffv - Fy + Gy*ptemp + y_py;
                        
      dV        = -CD\f;
      du        = dV(1:Nu);
      dv        = dV(Nu+1:Nu+Nv);
 
      V         = V + dV;
      uh        = uh + du;
      vh        = vh + dv;

%       error_nonlinear = norm(f)/norm(Om.*Vn/dt);
      error_nonlinear = norm(f,inf);
      
      i = i + 1;
 
    
end

if (i>= nonlinear_maxit || error_nonlinear>nonlinear_acc || isnan(error_nonlinear))
    error(['nonlinear iterations did not converge in ' num2str(nonlinear_maxit) ' iterations']);
end

R   = [uh; vh];

% evaluate divergence BC at endpoint, necessary for PC
t = tn + dt;
boundary_conditions;
interpolate_bc;
operator_bc_divergence;
t = tn;

% the factor of (1/2) for IM is taken into account by changing dt
% leading to second order pressure
dtn = dt;
dt  = 0.5*dt;
pressure_correction;
dt  = dtn;

p_old  = p;
p      = p + dp;


uh_old = uh;
vh_old = vh;

uh     = V(1:Nu);
vh     = V(Nu+1:Nu+Nv);


%% we expect a change in kinetic energy according to
% Vtemp = [utemp;vtemp];
% % 
% sum(Om.*(R-Vn).*(Vn+V))/2
% % 
% sum( (R-Vn).*(Om.*Vtemp-0.25*dt*(G*dp)))
% %
% dt*( - sum([yConv_ux + yConv_uy + Gx*pn; yConv_vx + yConv_vy + Gy*pn].*Vtemp) - sum(R.*(G*dp)/4))
% 

% utemp_p = Bup*(Au_ux * utemp + yAu_ux);
% vtemp_p = Bvp*(Av_vy * vtemp + yAv_vy);
% Vtemp_p = sqrt(utemp_p.^2 + vtemp_p.^2);

% 
% %
% dt*( - sum([yConv_ux + yConv_uy; yConv_vx + yConv_vy].*Vtemp) - sum(R.*(G*(pn+p))/4) )
% % ??? expect the minus for the first term to disappear here:
% 0.25*dt*( -sum((Vtemp_p.^2).*(M*R)) - sum(R.*(G*(pn+p))) )
% %
% 0.25*dt*( -sum((Vtemp_p.^2).*(M*R)) + sum((pn+p).*(M*R)) )
% %
% 0.25*dt*sum( (-Vtemp_p.^2 + (pn+p)).*(M*R) )
% %
% (1/8)*(dt^2)*sum(Om_inv.*( (G*(Vtemp_p.^2 - (p+pn))) .*(G*dp)) )
% %
% (1/8)*(dt^2)*sum(Om_inv.*( (G*(Vtemp_p.^2)).*(G*dp) + (G*pn).^2 - (G*p).^2) )

% k_error_linear(n) = (1/8)*(dt^2)*sum(Om_inv.*((G*pn).^2 - (G*p).^2) );
% k_error_nonlinear(n) = (1/8)*(dt^2)*sum(Om_inv.*( (G*(Vtemp_p.^2)).*(G*dp) ) );


% (dt^2/8)*(sum(Om_inv.*(G*p_old).^2) - sum(Om_inv.*(G*p).^2))
% max(M*(Vhn+V)/2)