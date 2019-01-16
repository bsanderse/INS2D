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


% coefficients for extrapolation of convective velocity
% Vandermonde matrix, in case of constant time step
% we could precompute the coefficients at the start of the time loop
VM_i    = method_startup_no-1:-1:0; 
VM1_row = c2+VM_i;
VM1     = vander(VM1_row)';
VM_rhs  = zeros(method_startup_no,1);
VM_rhs(end) = 1;

ep1     = VM1\VM_rhs;
 
% ep1     = [0;0;0;0];
% for non-constant time step:
% VM_i    = [0 dt(n-1) dt(n-2) dt(n-3)];
% VM1_row = c1*dt;

cV    = V_ep*ep1;
cu2   = cV(1:Nu);
cv2   = cV(Nu+1:Nu+Nv);


% initialization
f1     = zeros(Nu+Nv,1);
Z1     = spalloc(Np,Np,0);

% velocities at time level n
uhn    = uh;
vhn    = vh;
Vn     = V;
pn     = p;
tn     = t;

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


     
           
      %% rhs 1

      
      % convective terms, u-component
      uIux       = Iu_ux*cu2 + yIu_ux; % convective velocity, u_bar
      yConv_ux   = Cux*( uIux.*yAu_ux );

      vIuy       = Iv_uy*cv2 + yIv_uy; % convective velocity, v_bar
      yConv_uy   = Cuy*( vIuy.*yAu_uy );

      % convective terms, v-component
      uIvx       = Iu_vx*cu2 + yIu_vx; % convective velocity, u_bar  
      yConv_vx   = Cvx*( uIvx.*yAv_vx );

      vIvy       = Iv_vy*cv2 + yIv_vy; % convective velocity, v_bar
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

      
      % du/dt = k(u), u_i = u^n + dt*...
      ku2        = -yConv_ux - yConv_uy + yDiffu + Fx;
      kv2        = -yConv_vx - yConv_vy + yDiffv + Fy;
      
      % assemble rhs 1
      f2(1:Nu)        =  Omu.*uhn + dt*(a21*ku1 + a22*ku2 - y_px); 
      f2(Nu+1:Nu+Nv)  =  Omv.*vhn + dt*(a21*kv1 + a22*kv2 - y_py);
    

      %% assemble total right-hand side, this is -1*residual
      fM        =  -yM;

     
      %% total matrix   
      diag1 = [Omu;Omv];
      CD    = spdiags(diag1,0,Nu+Nv,Nu+Nv);
      CD    = CD + dt*a22*CD_1; 
            
    
      Z         = [CD dt*G; M Z1];        
      f         = [f2;fM];

      % solve
      dq        = Z\f;

      u2        = dq(1:Nu);
      v2        = dq(Nu+1:Nu+Nv);
      p2        = dq(Nu+Nv+1:Nu+Nv+Np);      


      
% CN is stiffly accurate, so:
uh = u2;
vh = v2;
V  = [uh;vh];

% following Hairer: second order (surprisingly?)
p  = (p2 - a21*pn)/a22;

% own method: first order
% p  = p2/c2;

% check energy-conservation of the system
% 0.5*sum(Omu.*uh.^2) + 0.5*sum(Omv.*vh.^2) - k(n-1)
% (dt^2/8)*( -sum(Om_inv.*([ku1;kv1] - G*pn).^2) + sum(Om_inv.*([ku2;kv2] - G*(p2/0.5-pn)).^2) )


% store new velocity field in V_ep
V_ep(:,1:method_startup_no-1) = V_ep(:,2:method_startup_no);
V_ep(:,method_startup_no)     = V;

pressure_additional_solve;