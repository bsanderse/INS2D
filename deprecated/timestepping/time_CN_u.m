% solution of saddlepoint system,
% nonlinear_Newton iteration

% solve for U_i

% reversible in time (inviscid flow)

% RK coefficients
a21    = 1/2;
a22    = 1/2;
b1     = 1/2;
b2     = 1/2;
c2     = 1;


f2     = zeros(Nu+Nv,1);

Z1     = spalloc(Np,Np,0);

if (n==2 || nonlinear_startingvalues==0)
    % initialization
    u2     = uh;
    v2     = vh;
    p2     = p;
else
    % better initialization using values from previous stage
    s        = 1;
    c_list   = [0 1];
    n_init   = length(c_list);
    id       = eye(n_init);
    for is=1:n_init
        poly_coeff     = polyfit(c_list,id(is,:),n_init-1);
        coeff(is,1:s)  = polyval(poly_coeff,1+c2);
    end
    u2_temp     = coeff(1,1)*uhn + coeff(2,1)*uh;
    v2_temp     = coeff(1,1)*vhn + coeff(2,1)*vh;
    p2_temp     = coeff(1,1)*pn + coeff(2,1)*p;
    u2 = u2_temp;
    v2 = v2_temp;
    p2 = p2_temp;
end


if (nonlinear_build_matrix==0)
      diag1 = [Omu;Omv];
      CD    = spdiags(diag1,0,Nu+Nv,Nu+Nv);
end

% velocities at time level n
uhn    = uh;
vhn    = vh;
Vn     = V;
pn     = p;
tn     = t;

error_nonlinear    = zeros(nonlinear_maxit,1);
error_nonlinear(1) = 1;
% iterate to remove linearization error; necessary to obtain full
% conservation and reversibility; if not removed still fully conservative
% but not reversible, and maybe not 2nd order
i = 0;

% F1 based on V1
ku1 = -Cux*( (Iu_ux*uhn+yIu_ux).*(Au_ux*uhn+yAu_ux) ) + ...
      -Cuy*( (Iv_uy*vhn+yIv_uy).*(Au_uy*uhn+yAu_uy) ) + ...
      Diffu*uhn + yDiffu;
kv1 = -Cvx*( (Iu_vx*uhn+yIu_vx).*(Av_vx*vhn+yAv_vx) ) + ...
      -Cvy*( (Iv_vy*vhn+yIv_vy).*(Av_vy*vhn+yAv_vy) ) + ...
      Diffv*vhn + yDiffv;

% evaluate BC and force at intermediate time
t = tn + c2*dt;
boundary_conditions;
interpolate_bc;
operator_bc_divergence;
operator_bc_momentum;
force;
t = tn;


while (i<nonlinear_maxit)
 
      i = i+1;
    
    
      % iterate until u^(n+1) = u^i = u^(i-1); this is different from
      % u^(n)!
    
           
      %% rhs 1

      
      % convective terms, u-component
      uIux       = Iu_ux*u2 + yIu_ux; % convective velocity, u_bar
      uAux       = Au_ux*u2 + yAu_ux;
      yConv_ux   = Cux*( uIux.*uAux );

      vIuy       = Iv_uy*v2 + yIv_uy; % convective velocity, v_bar
      uAuy       = Au_uy*u2 + yAu_uy;
      yConv_uy   = Cuy*( vIuy.*uAuy );

      % convective terms, v-component
      uIvx       = Iu_vx*u2 + yIu_vx; % convective velocity, u_bar  
      vAvx       = Av_vx*v2 + yAv_vx;
      yConv_vx   = Cvx*( uIvx.*vAvx );

      vIvy       = Iv_vy*v2 + yIv_vy; % convective velocity, v_bar
      vAvy       = Av_vy*v2 + yAv_vy;   
      yConv_vy   = Cvy*( vIvy.*vAvy );



      if (nonlinear_build_matrix==1 && (simplified_Newton==0 || i==1))


          % BLOCKS 11 - 12 - 21 - 22
          C1         = Cux*spdiags(uIux,0,N1,N1);   
          C2         = Cux*spdiags(uAux,0,N1,N1)*nonlinear_Newton;
          Conv_ux_11 = C1*Au_ux + C2*Iu_ux;

          C1         = Cuy*spdiags(vIuy,0,N2,N2);   
          C2         = Cuy*spdiags(uAuy,0,N2,N2)*nonlinear_Newton;
          Conv_uy_11 = C1*Au_uy;
          Conv_uy_12 = C2*Iv_uy;

          C1         = Cvx*spdiags(uIvx,0,N3,N3);
          C2         = Cvx*spdiags(vAvx,0,N3,N3)*nonlinear_Newton;
          Conv_vx_21 = C2*Iu_vx;
          Conv_vx_22 = C1*Av_vx;

          C1         = Cvy*spdiags(vIvy,0,N4,N4);
          C2         = Cvy*spdiags(vAvy,0,N4,N4)*nonlinear_Newton;       
          Conv_vy_22 = C1*Av_vy + C2*Iv_vy;


          % assemble matrix 1
          
          CD_11      = - Diffu + Conv_ux_11 + Conv_uy_11;
          CD_12      = spalloc(Nu,Nv,0) + Conv_uy_12;
          CD_21      = spalloc(Nv,Nu,0) + Conv_vx_21;
          CD_22      = - Diffv + Conv_vx_22 + Conv_vy_22; 

          CD_1       = [CD_11 CD_12; CD_21 CD_22];

      end

      
      % du/dt = k(u), u_i = u^n + dt*...
      ku2        = -yConv_ux - yConv_uy + Diffu*u2 + yDiffu + Fx;
      kv2        = -yConv_vx - yConv_vy + Diffv*v2 + yDiffv + Fy;
      
      % assemble rhs 1
      f2(1:Nu)        =  Omu.*(uhn-u2) + dt*(a21*ku1 + a22*ku2 - (Gx*p2 + y_px)); 
      f2(Nu+1:Nu+Nv)  =  Omv.*(vhn-v2) + dt*(a21*kv1 + a22*kv2 - (Gy*p2 + y_py));
    

      %% assemble total right-hand side, this is -1*residual
      fM        = -M*[u2;v2] - yM;
      f         = [f2;fM];

      error_nonlinear(i) = norm(f,inf);
      if (error_nonlinear(i) <= nonlinear_acc || ...
          error_nonlinear(i)/error_nonlinear(1)<= nonlinear_relacc)
          break;
      end
      
      %% total matrix       
      if (nonlinear_build_matrix==1 && (simplified_Newton==0 || i==1))
          diag1 = [Omu;Omv];
          CD    = spdiags(diag1,0,Nu+Nv,Nu+Nv);
          CD    = CD + dt*a22*CD_1; 
      end
            
      if (simplified_Newton==1)
          if (i==1)
              Z         = [CD dt*G; M Z1];      
%               [LZ,UZ]   = lu(Z);
          end
          % solve
%           dq        = UZ\(LZ\f);          % LU decomposition is too slow
          dq        = Z\f;   
      else
          Z         = [CD dt*G; M Z1];        
          dq        = Z\f;                  % this is generally faster than LU decomposition!
      end

      

      du2       = dq(1:Nu);
      dv2       = dq(Nu+1:Nu+Nv);
      dp2       = dq(Nu+Nv+1:Nu+Nv+Np);


      u2        = u2 + du2;
      v2        = v2 + dv2;
      p2        = p2 + dp2;
        
 
      % note that dp does not necessarily go to zero
      % divide the residual by a norm based on rhs ?
%       k_norm = [Om.*Vn;Om.*Vn];
%       error_nonlinear = norm(f,inf) /norm(k_norm,inf);
%       error_nonlinear = max(norm(f),norm(dq))/norm(k_norm)
      
end


if ((error_nonlinear(i)>nonlinear_acc && error_nonlinear(i)/error_nonlinear(1)> nonlinear_relacc) ...
     || i>= nonlinear_maxit || isnan(error_nonlinear(i)))
    error(['nonlinear iterations did not converge in ' num2str(nonlinear_maxit) ' iterations']);
    error_nonlinear(i)
    error_nonlinear(i)/error_nonlinear(1)
end
nonlinear_its(n) = i-1;


% CN is stiffly accurate, so:
uh = u2;
vh = v2;
V  = [uh;vh];

% following Hairer: second order (surprisingly?)
p  = (p2 - a21*pn)/a22;

% own method: first order
% p  = p2/c2;



% check energy-conservation of the system
% we need to use the 'pressure' p, not the lagrange multiplier p2
% calculated in the nonlinear routine
% 0.5*sum(Omu.*uh.^2) + 0.5*sum(Omv.*vh.^2) - k(n-1)
% (dt^2/8)*( -sum(Om_inv.*([ku1;kv1] - G*pn).^2) + sum(Om_inv.*([ku2;kv2] - G*(p2/0.5-pn)).^2) )
% keyboard

% store velocity fields if this method is used as startup 
% for linear extrapolation
if (method_temp == 192)
    V_ep(:,n) = V;  
end

pressure_additional_solve;