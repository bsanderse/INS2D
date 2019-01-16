% fully (energy) conservative time stepping
% solution of saddlepoint system,
% nonlinear_Newton iteration

% solve for U_i

% reversible in time (inviscid flow)

% RK coefficients
a11    = 1/2;
b1     = 1;
c1     = 1/2;

a_inv_11 = 2;


f1     = zeros(Nu+Nv,1);

Z1     = spalloc(Np,Np,0);

if (n==2 || nonlinear_startingvalues==0)
    % initialization
    u1     = uh;
    v1     = vh;
    p1     = p;
else
    % better initialization using values from previous stage
    s        = 1;
    c_list   = [0 c1 1];
    n_init   = length(c_list);
    id       = eye(n_init);
    for is=1:n_init
        poly_coeff     = polyfit(c_list,id(is,:),n_init-1);
        coeff(is,1:s)  = polyval(poly_coeff,1+c1);
    end
    u1_temp     = coeff(1,1)*uhn + coeff(2,1)*u1 + coeff(3,1)*uh;
    v1_temp     = coeff(1,1)*vhn + coeff(2,1)*v1 + coeff(3,1)*vh;
    p1_temp     = coeff(1,1)*pn + coeff(2,1)*p1 + coeff(3,1)*p;
    u1 = u1_temp;
    v1 = v1_temp;
    p1 = p1_temp;
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


% evaluate BC and force at intermediate times 
t = tn + c1*dt;
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
      uIux       = Iu_ux*u1 + yIu_ux; % convective velocity, u_bar
      uAux       = Au_ux*u1 + yAu_ux;
      yConv_ux   = Cux*( uIux.*uAux );

      vIuy       = Iv_uy*v1 + yIv_uy; % convective velocity, v_bar
      uAuy       = Au_uy*u1 + yAu_uy;
      yConv_uy   = Cuy*( vIuy.*uAuy );

      % convective terms, v-component
      uIvx       = Iu_vx*u1 + yIu_vx; % convective velocity, u_bar  
      vAvx       = Av_vx*v1 + yAv_vx;
      yConv_vx   = Cvx*( uIvx.*vAvx );

      vIvy       = Iv_vy*v1 + yIv_vy; % convective velocity, v_bar
      vAvy       = Av_vy*v1 + yAv_vy;   
      yConv_vy   = Cvy*( vIvy.*vAvy );

      if (order4==1)
          
          uIux3      = Iu_ux3*u1 + yIu_ux3;
          uAux3      = Au_ux3*u1 + yAu_ux3;   
          yConv_ux_3 = Cux3*( uIux3.*uAux3 );

          vIuy3      = Iv_uy3*v1 + yIv_uy3;
          uAuy3      = Au_uy3*u1 + yAu_uy3;
          yConv_uy_3 = Cuy3*( vIuy3.*uAuy3 );
          
          uIvx3      = Iu_vx3*u1 + yIu_vx3;
          vAvx3      = Av_vx3*v1 + yAv_vx3;
          yConv_vx_3 = Cvx3*( uIvx3.*vAvx3 );

          vIvy3      = Iv_vy3*v1 + yIv_vy3;
          vAvy3      = Av_vy3*v1 + yAv_vy3;
          yConv_vy_3 = Cvy3*( vIvy3.*vAvy3 );
          
      end
      

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


          if (order4==1)
              C1         = Cux3*spdiags(uIux3,0,length(uIux3),length(uIux3));   
              C2         = Cux3*spdiags(uAux3,0,length(uAux3),length(uAux3))*nonlinear_Newton;
              Conv_ux_11_3 = C1*Au_ux3 + C2*Iu_ux3;
              
              C1         = Cuy3*spdiags(vIuy3,0,length(vIuy3),length(vIuy3));   
              C2         = Cuy3*spdiags(uAuy3,0,length(vIuy3),length(vIuy3))*nonlinear_Newton;
              Conv_uy_11_3 = C1*Au_uy3;
              Conv_uy_12_3 = C2*Iv_uy3;
              
              C1         = Cvx3*spdiags(uIvx3,0,length(uIvx3),length(uIvx3));
              C2         = Cvx3*spdiags(vAvx3,0,length(vAvx3),length(vAvx3))*nonlinear_Newton;
              Conv_vx_21_3 = C2*Iu_vx3;
              Conv_vx_22_3 = C1*Av_vx3;
              
              C1         = Cvy3*spdiags(vIvy3,0,length(vIvy3),length(vIvy3));
              C2         = Cvy3*spdiags(vAvy3,0,length(vAvy3),length(vAvy3))*nonlinear_Newton;       
              Conv_vy_22_3 = C1*Av_vy3 + C2*Iv_vy3;     
          end
          
          % assemble matrix 1
          
          if (order4==0)
              CD_11  = - Diffu + Conv_ux_11 + Conv_uy_11;
              CD_12  = spalloc(Nu,Nv,0) + Conv_uy_12;
              CD_21  = spalloc(Nv,Nu,0) + Conv_vx_21;
              CD_22  = - Diffv + Conv_vx_22 + Conv_vy_22;
          elseif (order4==1)
              CD_11  = - Diffu + alfa*Conv_ux_11 - Conv_ux_11_3 + alfa*Conv_uy_11 - Conv_uy_11_3;
              CD_12  = alfa*Conv_uy_12 - Conv_uy_12_3;
              CD_21  = alfa*Conv_vx_21 - Conv_vx_21_3;
              CD_22  = - Diffv + alfa*Conv_vx_22 - Conv_vx_22_3 + alfa*Conv_vy_22 - Conv_vy_22_3;  
          end
              
          CD_1       = [CD_11 CD_12; CD_21 CD_22];
%           max2d(abs(CD_11+CD_11'))
%           max2d(abs(CD_22+CD_22'))
      end

      
      % du/dt = k(u), u_i = u^n + dt*a11*ku1
      if (order4==0)
          ku1        = -yConv_ux - yConv_uy + Diffu*u1 + yDiffu + Fx;% - Gx*p1 - y_px_1;
          kv1        = -yConv_vx - yConv_vy + Diffv*v1 + yDiffv + Fy;% - Gy*p1 - y_py_1;      
      elseif (order4==1)
          ku1        = - (alfa*yConv_ux - yConv_ux_3) - (alfa*yConv_uy - yConv_uy_3) + ...
                       Diffu*u1 + yDiffu + Fx;
          kv1        = - (alfa*yConv_vx - yConv_vx_3) - (alfa*yConv_vy - yConv_vy_3) + ...
                       Diffv*v1 + yDiffv + Fy;
      end
      % assemble rhs 1
      f1(1:Nu)        =  Omu.*(uhn-u1) + dt*(a11*ku1  - c1*(Gx*p1 + y_px)); 
      f1(Nu+1:Nu+Nv)  =  Omv.*(vhn-v1) + dt*(a11*kv1  - c1*(Gy*p1 + y_py));
    

      %% assemble total right-hand side, this is -1*residual
      fM        = -M*[u1;v1] - yM;
      f         = [f1;fM];

      error_nonlinear(i) = norm(f,inf);
      if (error_nonlinear(i) <= nonlinear_acc || ...
          error_nonlinear(i)/error_nonlinear(1)<= nonlinear_relacc)
          break;
      end
      
      %% total matrix       
      if (nonlinear_build_matrix==1 && (simplified_Newton==0 || i==1))
          diag1 = [Omu;Omv];
          CD    = spdiags(diag1,0,Nu+Nv,Nu+Nv);
          CD    = CD + dt*a11*CD_1; 
      end
            
      if (simplified_Newton==1)
          if (i==1)
              Z         = [CD c1*dt*G; M Z1];      
%               [LZ,UZ]   = lu(Z);
          end
          % solve
%           dq        = UZ\(LZ\f);          % LU decomposition is too slow
          dq        = Z\f;   
      else
          Z         = [CD c1*dt*G; M Z1];        
          dq        = Z\f;                  % this is generally faster than LU decomposition!
      end

      

      du1       = dq(1:Nu);
      dv1       = dq(Nu+1:Nu+Nv);
      dp1       = dq(Nu+Nv+1:Nu+Nv+Np);


      u1        = u1 + du1;
      v1        = v1 + dv1;
      p1        = p1 + dp1;
        
 
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


ku1_d = Diffu*u1 + yDiffu;
kv1_d = Diffv*v1 + yDiffv;
if (order4==0)
    ku1_c = -Cux*( (Iu_ux*u1+yIu_ux).*(Au_ux*u1+yAu_ux) ) + ...
            -Cuy*( (Iv_uy*v1+yIv_uy).*(Au_uy*u1+yAu_uy) );
    kv1_c = -Cvx*( (Iu_vx*u1+yIu_vx).*(Av_vx*v1+yAv_vx) ) + ...
            -Cvy*( (Iv_vy*v1+yIv_vy).*(Av_vy*v1+yAv_vy) );
elseif (order4==1)
    ku1_c = - alfa*Cux*(  (Iu_ux*u1+yIu_ux).*  (Au_ux*u1+yAu_ux) ) + ...
                 + Cux3*( (Iu_ux3*u1+yIu_ux3).*(Au_ux3*u1+yAu_ux3) ) ...
            - alfa*Cuy*(  (Iv_uy*v1+yIv_uy).*  (Au_uy*u1+yAu_uy) ) + ...
                 + Cuy3*( (Iv_uy3*v1+yIv_uy3).*(Au_uy3*u1+yAu_uy3) );
    kv1_c = - alfa*Cvx*(  (Iu_vx*u1+yIu_vx).*  (Av_vx*v1+yAv_vx) ) + ...
                 + Cvx3*( (Iu_vx3*u1+yIu_vx3).*(Av_vx3*v1+yAv_vx3) ) ...
            - alfa*Cvy*(  (Iv_vy*v1+yIv_vy).*  (Av_vy*v1+yAv_vy) ) + ...
                 + Cvy3*( (Iv_vy3*v1+yIv_vy3).*(Av_vy3*v1+yAv_vy3) );
    
end
uh    = uhn + dt*Omu_inv.*( b1*( ku1_d + ku1_c + Fx - y_px) - b1*a_inv_11*c1*Gx*p1);
vh    = vhn + dt*Omv_inv.*( b1*( kv1_d + kv1_c + Fy - y_py) - b1*a_inv_11*c1*Gy*p1);
 

% update with U1 and U2 
% with the following formulation we don't have to recompute the right hand side with
% the new velocity vector, but simply use the inverse of the Butcher
% tableau; however, this is more sensitive to round-off error accumulation
% uh = uhn + b1*a_inv_11*(u1-uhn);
% vh = vhn + b1*a_inv_11*(v1-vhn);
V  = [uh;vh];


if (BC_unsteady == 1)
    % make V satisfy the incompressibility constraint at n+1
    t = tn + dt;
    boundary_conditions;
    interpolate_bc;
    operator_bc_divergence;
    t = tn;

    f       = (1/dt)*(M*V + yM);

    pressure_poisson;

    V       = V - dt*Om_inv.*(G*dp);
    uh      = V(1:Nu);
    vh      = V(Nu+1:Nu+Nv);

%     p       = pn + dp;
    
end

% max(abs(M*V+yM))

% if (max(abs(M*V+yM))>nonlinear_acc)
%     keyboard;
% end 
   
% psi1  = a_inv_11*c1*p1;
% p     = psi1;
p = p1;
% p = 2*p1-pn;

% store velocity fields if this method is used as startup 
% for linear extrapolation
if (method_temp == 62)
    V_ep(:,n) = V;  
end

pressure_additional_solve;