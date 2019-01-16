% fully (energy) conservative time stepping
% solution of saddlepoint system,
% nonlinear_Newton iteration

% reversible in time (inviscid flow)

% RK coefficients
a11    = 1/2;
b1     = 1;
c1     = 1/2;

a_inv_11 = 2;

f      = zeros(Nu+Nv,1);

% velocities at time level n
uhn    = uh;
vhn    = vh;
Vn     = V;
pn     = p;
tn     = t;


% initialization
u1     = uh;
v1     = vh;
ku1    = zeros(Nu,1);
kv1    = zeros(Nv,1);
p1     = p;


if (nonlinear_build_matrix==0)
      diag1 = [Omu;Omv];
      CD    = spdiags(diag1,0,Nu+Nv,Nu+Nv);
end

error_nonlinear = 1;
% iterate to remove linearization error; necessary to obtain full
% conservation and reversibility; if not removed still fully conservative
% but not reversible, and maybe not 2nd order
i=1;


% evaluate BC and force at intermediate time 
t = tn + c1*dt;
if (BC_unsteady == 1)
    boundary_conditions;
    interpolate_bc;
    operator_bc_divergence;
    operator_bc_momentum;        
end
force;
t = tn;


while (error_nonlinear > nonlinear_acc  && i<nonlinear_maxit)
 
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

      % assemble rhs
      f(1:Nu)   =  Omu.*ku1 + yConv_ux + yConv_uy + ...
                   -Diffu*u1 - yDiffu - Fx + Gx*pn + y_px;
                        
      f(Nu+1:Nu+Nv)  =  Omv.*kv1 + yConv_vx + yConv_vy + ...
                        -Diffv*v1 - yDiffv - Fy + Gy*pn + y_py;

      if (nonlinear_build_matrix==1)


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


          % assemble matrix
          
          CD_11      = - Diffu + Conv_ux_11 + Conv_uy_11;
          CD_12      = spalloc(Nu,Nv,0) + Conv_uy_12; % in this way we can comment out convection
          CD_21      = spalloc(Nv,Nu,0) + Conv_vx_21;
          CD_22      = - Diffv  + Conv_vx_22 + Conv_vy_22; 

          CD_1       = [CD_11 CD_12; CD_21 CD_22];

          diag1      = [Omu;Omv];
          CD         = spdiags(diag1,0,Nu+Nv,Nu+Nv);
          CD         = CD + dt*a11*CD_1; 
          
      end


      %% assemble total right-hand side, this is -1*residual

      % solve
      dq        = -CD\f;
 
      dku1      = dq(1:Nu);
      dkv1      = dq(Nu+1:Nu+Nv);

      ku1       = ku1 + dku1;
      kv1       = kv1 + dkv1;
      
      u1        = uhn + dt*a11*ku1;
      v1        = vhn + dt*a11*kv1;
     
      
      % note that dp does not necessarily go to zero
      % divide the residual by a norm based on rhs ?
      k_norm = [Om.*Vn;Om.*Vn];
      error_nonlinear = norm(f,inf);%/norm(k_norm,inf);
%       error_nonlinear = norm(f,inf);
%       error_nonlinear = max(norm(f),norm(dq))/norm(k_norm)
      i = i+1;
%       pause;
      
end

if (i>= nonlinear_maxit || error_nonlinear>nonlinear_acc || isnan(error_nonlinear))
    error(['nonlinear iterations did not converge in ' num2str(nonlinear_maxit) ' iterations']);
end

uh      = uhn + dt*b1*ku1;
vh      = vhn + dt*b1*kv1;
V       = [uh;vh];

% make V satisfy the incompressibility constraint at n+1 with a
% projection step
if (BC_unsteady==1)
    t = tn + dt;
    boundary_conditions;
    interpolate_bc;
    operator_bc_divergence;
    t = tn;
end

f       = M*V + yM;

pressure_poisson;
p       = pn + dp;

V       = V - Om_inv.*(G*dp);
uh      = V(1:Nu);
vh      = V(Nu+1:Nu+Nv);


pressure_additional_solve;