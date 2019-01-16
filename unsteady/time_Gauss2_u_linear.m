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

% coefficients for extrapolation of convective velocity
% Vandermonde matrix, in case of constant time step
% we could precompute the coefficients at the start of the time loop
VM_i    = method_startup_no-1:-1:0; 
VM1_row = c1+VM_i;
VM1     = vander(VM1_row)';
VM_rhs  = zeros(method_startup_no,1);
VM_rhs(end) = 1;

ep1     = VM1\VM_rhs;
 
% ep1     = [0;0;0;0];
% for non-constant time step:
% VM_i    = [0 dt(n-1) dt(n-2) dt(n-3)];
% VM1_row = c1*dt;

cV1    = V_ep*ep1;
cu     = cV1(1:Nu);
cv     = cV1(Nu+1:Nu+Nv);


% initialization
f1     = zeros(Nu+Nv,1);
Z1     = spalloc(Np,Np,0);

% velocities at time level n
uhn    = uh;
vhn    = vh;
Vn     = V;
pn     = p;
tn     = t;


% evaluate BC and force at intermediate times 
t = tn + c1*dt;
boundary_conditions;
interpolate_bc;
operator_bc_divergence;
operator_bc_momentum;
force;
t = tn;
    
           
      %% rhs 1

      
      % convective terms, u-component
      uIux       = Iu_ux*cu + yIu_ux; % convective velocity, u_bar
      yConv_ux   = Cux*( uIux.*yAu_ux );

      vIuy       = Iv_uy*cv + yIv_uy; % convective velocity, v_bar
      yConv_uy   = Cuy*( vIuy.*yAu_uy );

      % convective terms, v-component
      uIvx       = Iu_vx*cu + yIu_vx; % convective velocity, u_bar  
      yConv_vx   = Cvx*( uIvx.*yAv_vx );

      vIvy       = Iv_vy*cv + yIv_vy; % convective velocity, v_bar
      yConv_vy   = Cvy*( vIvy.*yAv_vy );


      % BLOCKS 11 - 12 - 21 - 22
      C1         = Cux*spdiags(uIux,0,N1,N1);   
      Conv_ux_11 = C1*Au_ux;

      C1         = Cuy*spdiags(vIuy,0,N2,N2);   
      Conv_uy_11 = C1*Au_uy;

      C1         = Cvx*spdiags(uIvx,0,N3,N3);
      Conv_vx_22 = C1*Av_vx;

      C1         = Cvy*spdiags(vIvy,0,N4,N4);
      Conv_vy_22 = C1*Av_vy;


      % assemble matrix 1

      CD_11      = - Diffu + Conv_ux_11 + Conv_uy_11;
      CD_12      = spalloc(Nu,Nv,0);
      CD_21      = spalloc(Nv,Nu,0);
      CD_22      = - Diffv + Conv_vx_22 + Conv_vy_22; 

      CD_1       = [CD_11 CD_12; CD_21 CD_22];


      
      % du/dt = k(u), u_i = u^n + dt*a11*ku1
      ku1        = -yConv_ux - yConv_uy + yDiffu + Fx;% - Gx*p1 - y_px_1;
      kv1        = -yConv_vx - yConv_vy + yDiffv + Fy;% - Gy*p1 - y_py_1;      
      
      % assemble rhs 1
      f1(1:Nu)        =  Omu.*uhn + dt*(a11*ku1  - c1*y_px); 
      f1(Nu+1:Nu+Nv)  =  Omv.*vhn + dt*(a11*kv1  - c1*y_py);
    

      %% assemble total right-hand side, this is -1*residual
      fM         =  -yM;

      
      %% total matrix       
      diag1 = [Omu;Omv];
      CD    = spdiags(diag1,0,Nu+Nv,Nu+Nv);
      CD    = CD + dt*a11*CD_1; 
            
    
      Z         = [CD c1*dt*G; M Z1];        
      f         = [f1;fM];

      % solve
      dq        = Z\f;

      u1        = dq(1:Nu);
      v1        = dq(Nu+1:Nu+Nv);
      p1        = dq(Nu+Nv+1:Nu+Nv+Np);
  

ku1_d = Diffu*u1 + yDiffu;
kv1_d = Diffv*v1 + yDiffv;
ku1_c = -Cux*( (Iu_ux*cu+yIu_ux).*(Au_ux*u1+yAu_ux) ) + ...
        -Cuy*( (Iv_uy*cv+yIv_uy).*(Au_uy*u1+yAu_uy) );
kv1_c = -Cvx*( (Iu_vx*cu+yIu_vx).*(Av_vx*v1+yAv_vx) ) + ...
        -Cvy*( (Iv_vy*cv+yIv_vy).*(Av_vy*v1+yAv_vy) );
uh    = uhn + dt*Omu_inv.*( b1*( ku1_d + ku1_c + Fx - y_px) - b1*a_inv_11*c1*Gx*p1);
vh    = vhn + dt*Omv_inv.*( b1*( kv1_d + kv1_c + Fy - y_py) - b1*a_inv_11*c1*Gy*p1);
       
% or update with U1 and U2 
% with the following formulation we don't have to recompute the right hand side with
% the new velocity vector, but simply use the inverse of the Butcher
% tableau
% uh = uhn + b1*a_inv_11*(u1-uhn);
% vh = vhn + b1*a_inv_11*(v1-vhn);
V  = [uh;vh];

% p  = pn + b1*a_inv_11*(p1-pn);


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

    p       = pn + dp;
    
end

% psi1  = a_inv_11*c1*p1;
% p     = psi1;
p = p1;

% store new velocity field in V_ep
V_ep(:,1:method_startup_no-1) = V_ep(:,2:method_startup_no);
V_ep(:,method_startup_no)     = V;

pressure_additional_solve;