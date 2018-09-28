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

if (n==2)
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


grad = c1*dt*G;

if (nonlinear_build_matrix==0 && n==2)
      diag1 = [Omu;Omv];
      CD    = spdiags(diag1,0,Nu+Nv,Nu+Nv);
      CD    = diag(1./diag(CD));     
      % poisson-like matrix AS
      AS    = M*CD*grad;
      [LS,US] = lu(AS); 
end

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


% evaluate BC and force at intermediate times 
t = tn + c1*dt;
boundary_conditions;
interpolate_bc;
operator_bc_divergence;
operator_bc_momentum;
force;
t = tn;




while (error_nonlinear > nonlinear_acc && i<nonlinear_maxit)
 
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



      if (nonlinear_build_matrix==1)


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

      
      % du/dt = k(u), u_i = u^n + dt*a11*ku1
      ku1        = -yConv_ux - yConv_uy + Diffu*u1 + yDiffu + Fx;% - Gx*p1 - y_px_1;
      kv1        = -yConv_vx - yConv_vy + Diffv*v1 + yDiffv + Fy;% - Gy*p1 - y_py_1;      
      
      % assemble rhs 1
      f1(1:Nu)        =  Omu.*(uhn-u1) + dt*(a11*ku1  - c1*(Gx*p1 + y_px)); 
      f1(Nu+1:Nu+Nv)  =  Omv.*(vhn-v1) + dt*(a11*kv1  - c1*(Gy*p1 + y_py));
    

      %% assemble total right-hand side, this is -1*residual
      fM        = -M*[u1;v1] - yM;

      
      %% total matrix       
      if (nonlinear_build_matrix==1)
          diag1 = [Omu;Omv];
          CD    = spdiags(diag1,0,Nu+Nv,Nu+Nv);
          CD    = CD + dt*a11*CD_1;        
          CD    = diag(1./diag(CD));
      end 
                
      % solve for velocity
      dq        = CD*f1;
     
      % poisson like for the pressure
      if (nonlinear_build_matrix==1)     
          AS        = M*CD*grad;
          dp        = AS\(M*dq - fM);
      else
          dp        = US\(LS\(M*dq-fM));
      end
      
      % update velocity with pressure change
      dq        = dq - CD*grad*dp;

      du1       = dq(1:Nu);
      dv1       = dq(Nu+1:Nu+Nv);
      u1        = u1 + du1;
      v1        = v1 + dv1;
      p1        = p1 + dp;
        
      f         = [f1;fM];
 
      % note that dp does not necessarily go to zero
      % divide the residual by a norm based on rhs ?
      k_norm = [Om.*Vn;Om.*Vn];
%       error_nonlinear = norm(f,inf) /norm(k_norm,inf);
      error_nonlinear = norm(f,inf);
%       error_nonlinear = max(norm(f),norm(dq))/norm(k_norm)
      i = i+1;
      
end


if (i>= nonlinear_maxit || error_nonlinear>nonlinear_acc || isnan(error_nonlinear))
    error(['nonlinear iterations did not converge in ' num2str(nonlinear_maxit) ' iterations']);
end
nonlinear_its(n) = i-1;


% update with F1 and F2 
% uh   = uhn + dt*Omu_inv.*(b1*ku1 + b2*ku2);
% vh   = vhn + dt*Omv_inv.*(b1*kv1 + b2*kv2);
% V    = [uh;vh];

% or update with U1 and U2 
% with the following formulation we don't have to recompute the right hand side with
% the new velocity vector, but simply use the inverse of the Butcher
% tableau
uh = uhn + b1*a_inv_11*(u1-uhn);
vh = vhn + b1*a_inv_11*(v1-vhn);
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

psi1  = a_inv_11*c1*p1;
p     = psi1;

% store velocity fields if this method is used as startup 
% for linear extrapolation
if (method_temp == 62)
    V_ep(:,n) = V;  
end

pressure_additional_solve;