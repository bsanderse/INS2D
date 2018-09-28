% Radau IIA method of order 3
% solution of saddlepoint system,
% nonlinear_Newton iteration

% solve for U_i

radau_type = type_list(jj); % 0: Radau IIA, 1: Radau IIB, 2: Radau IIAB

% RK coefficients for convection
if (radau_type==0)
    a11c    = 5/12;
    a12c    = -1/12;
    a21c    = 3/4;
    a22c    = 1/4;
else
    a11c    = 3/8;
    a12c    = -1/24;
    a21c    = 7/8;
    a22c    = 1/8;
end

% RK coefficients for diffusion
if (radau_type==1)
    a11d    = 3/8;
    a12d    = -1/24;
    a21d    = 7/8;
    a22d    = 1/8;
else
    a11d    = 5/12;
    a12d    = -1/12;
    a21d    = 3/4;
    a22d    = 1/4;
end

% b and c-coefficients are same for convection and diffusion
b1     = 3/4;
b2     = 1/4;

c1     = 1/3;
c2     = 1;

% inverse of A-coefficients
if (radau_type==1) 
    % convection:
    a_inv_11 = 3/2;
    a_inv_12 = 1/2;
    a_inv_21 = -21/2;
    a_inv_22 = 9/2;
else
    % diffusion:
    a_inv_11 = 3/2;
    a_inv_12 = 1/2;
    a_inv_21 = -9/2;
    a_inv_22 = 5/2;
end

% coefficients for extrapolation of convective velocity
% Vandermonde matrix, in case of constant time step
% we could precompute the coefficients at the start of the time loop
VM_i    = method_startup_no-1:-1:0; 
VM1_row = c1+VM_i;
VM2_row = c2+VM_i;
VM1     = vander(VM1_row)';
VM2     = vander(VM2_row)';
VM_rhs  = zeros(method_startup_no,1);
VM_rhs(end) = 1;

ep1     = VM1\VM_rhs;
ep2     = VM2\VM_rhs;

cV1    = V_ep*ep1;
cV2    = V_ep*ep2;
cu1    = cV1(1:Nu);
cv1    = cV1(Nu+1:Nu+Nv);
cu2    = cV2(1:Nu);
cv2    = cV2(Nu+1:Nu+Nv);


% initialization
f1     = zeros(Nu+Nv,1);
f2     = zeros(Nu+Nv,1);
Z1     = spalloc(Np,Nu+Nv,0);
Z2     = spalloc(2*Np,2*Np,0);


% update velocities at time level n
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

yIu_ux_1   = yIu_ux;
yAu_ux_1   = yAu_ux;
yIv_uy_1   = yIv_uy;
yAu_uy_1   = yAu_uy;
yIu_vx_1   = yIu_vx;
yAv_vx_1   = yAv_vx;
yIv_vy_1   = yIv_vy;
yAv_vy_1   = yAv_vy;
yDiffu_1   = yDiffu;
yDiffv_1   = yDiffv;
Fx_1       = Fx;
Fy_1       = Fy;
y_px_1     = y_px;
y_py_1     = y_py;
yM1        = yM;


t = tn + c2*dt;
boundary_conditions;
interpolate_bc;
operator_bc_divergence;
operator_bc_momentum;
force;
t = tn;

yIu_ux_2   = yIu_ux;
yAu_ux_2   = yAu_ux;
yIv_uy_2   = yIv_uy;
yAu_uy_2   = yAu_uy;
yIu_vx_2   = yIu_vx;
yAv_vx_2   = yAv_vx;
yIv_vy_2   = yIv_vy;
yAv_vy_2   = yAv_vy;
yDiffu_2   = yDiffu;
yDiffv_2   = yDiffv;
Fx_2       = Fx;
Fy_2       = Fy;
y_px_2     = y_px;
y_py_2     = y_py;
yM2        = yM;   

 
           
  %% rhs 1

      
  % convective terms, u-component
  uIux       = Iu_ux*cu1 + yIu_ux_1; % convective velocity, u_bar
  yConv_ux_1 = Cux*( uIux.*yAu_ux_1 );

  vIuy       = Iv_uy*cv1 + yIv_uy_1; % convective velocity, v_bar
  yConv_uy_1 = Cuy*( vIuy.*yAu_uy_1 );

  % convective terms, v-component
  uIvx       = Iu_vx*cu1 + yIu_vx_1; % convective velocity, u_bar  
  yConv_vx_1 = Cvx*( uIvx.*yAv_vx_1 );

  vIvy       = Iv_vy*cv1 + yIv_vy_1; % convective velocity, v_bar
  yConv_vy_1 = Cvy*( vIvy.*yAv_vy_1 );

  
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

    C_11      = Conv_ux_11 + Conv_uy_11;
    C_12      = spalloc(Nu,Nv,0);
    C_21      = spalloc(Nv,Nu,0);
    C_22      = Conv_vx_22 + Conv_vy_22; 

    C_1       = [C_11 C_12; C_21 C_22];
    D_1       = [Diffu spalloc(Nu,Nv,0); spalloc(Nv,Nu,0) Diffv];

      
  %% rhs 2
      
  % convective terms, u-component
  uIux       = Iu_ux*cu2 + yIu_ux_2; % convective velocity, u_bar
  yConv_ux_2 = Cux*( uIux.*yAu_ux_2 );

  vIuy       = Iv_uy*cv2 + yIv_uy_2; % convective velocity, v_bar
  yConv_uy_2 = Cuy*( vIuy.*yAu_uy_2 );

  % convective terms, v-component
  uIvx       = Iu_vx*cu2 + yIu_vx_2; % convective velocity, u_bar  
  yConv_vx_2 = Cvx*( uIvx.*yAv_vx_2 );

  vIvy       = Iv_vy*cv2 + yIv_vy_2; % convective velocity, v_bar
  yConv_vy_2 = Cvy*( vIvy.*yAv_vy_2 );     

   
  C1         = Cux*spdiags(uIux,0,N1,N1);   
  Conv_ux_33 = C1*Au_ux;

  C1         = Cuy*spdiags(vIuy,0,N2,N2);   
  Conv_uy_33 = C1*Au_uy;

  C1         = Cvx*spdiags(uIvx,0,N3,N3);
  Conv_vx_44 = C1*Av_vx;

  C1         = Cvy*spdiags(vIvy,0,N4,N4);
  Conv_vy_44 = C1*Av_vy;   

  % assemble matrix 2   
  C_33      = Conv_ux_33 + Conv_uy_33;
  C_34      = spalloc(Nu,Nv,0);
  C_43      = spalloc(Nv,Nu,0);
  C_44      = Conv_vx_44 + Conv_vy_44; 

  C_2       = [C_33 C_34; C_43 C_44];
  D_2       = [Diffu spalloc(Nu,Nv,0); spalloc(Nv,Nu,0) Diffv];
          

  % du/dt = k(u), u_i = u^n + dt*(a11*ku1 + a12*ku2)
  ku1_c        = -yConv_ux_1 - yConv_uy_1;% - Gx*p1 - y_px_1;
  kv1_c        = -yConv_vx_1 - yConv_vy_1;% - Gy*p1 - y_py_1;
  ku1_d        = yDiffu_1 + Fx_1;
  kv1_d        = yDiffv_1 + Fy_1;

  ku2_c        = -yConv_ux_2 - yConv_uy_2;% - Gx*p2 - y_px_2;
  kv2_c        = -yConv_vx_2 - yConv_vy_2;% - Gy*p2 - y_py_2;
  ku2_d        = yDiffu_2 + Fx_2;
  kv2_d        = yDiffv_2 + Fy_2;
      
      
  % assemble rhs 1
  f1(1:Nu)        =  Omu.*uhn + dt*(a11c*ku1_c + a12c*ku2_c + ...
                                    a11d*ku1_d + a12d*ku2_d + ...
                                     - c1*y_px_1);
  f1(Nu+1:Nu+Nv)  =  Omv.*vhn + dt*(a11c*kv1_c + a12c*kv2_c + ...
                                    a11d*kv1_d + a12d*kv2_d + ...
                                     - c1*y_py_1);

  % assemble rhs 2 
  f2(1:Nu)        =  Omu.*uhn + dt*(a21c*ku1_c + a22c*ku2_c + ...
                                    a21d*ku1_d + a22d*ku2_d + ...
                                     - c2*y_px_2);                        
  f2(Nu+1:Nu+Nv)  =  Omv.*vhn + dt*(a21c*kv1_c + a22c*kv2_c + ...
                                    a21d*kv1_d + a22d*kv2_d + ...
                                     - c2*y_py_2);

  %% assemble total right-hand side, this is -1*residual
  fM1       = - yM1;
  fM2       = - yM2;

      
  %% total matrix       
  diag1 = [Omu;Omv;Omu;Omv];
  CD    = spdiags(diag1,0,2*(Nu+Nv),2*(Nu+Nv));
  CD    = CD + dt*[a11c*C_1 - a11d*D_1  a12c*C_2 - a12d*D_2; ...
                   a21c*C_1 - a21d*D_1  a22c*C_2 - a22d*D_2]; 

            
  grad    = [c1*dt*G Z1'; Z1' c2*dt*G];     
  div     = [M Z1; Z1 M];        
  Z       = [CD grad; div Z2];        
  f       = [f1;f2;fM1;fM2];

  % solve
  dq      = Z\f;

  u1      = dq(1:Nu);
  v1      = dq(Nu+1:Nu+Nv);
  u2      = dq(Nu+Nv+1:2*Nu+Nv);
  v2      = dq(2*Nu+Nv+1:2*(Nu+Nv));
  p1      = dq(2*(Nu+Nv)+1:2*(Nu+Nv)+Np);
  p2      = dq(2*(Nu+Nv)+Np+1:2*(Nu+Nv)+2*Np);


if (radau_type==1 || radau_type==2)
    % we have to compute F1 and F2 based on V1 and V2
    ku1_d = Diffu*u1 + yDiffu_1;
    kv1_d = Diffv*v1 + yDiffv_1;
    ku2_d = Diffu*u2 + yDiffu_2;
    kv2_d = Diffv*v2 + yDiffv_2;
    ku1_c = -Cux*( (Iu_ux*cu1+yIu_ux_1).*(Au_ux*u1+yAu_ux_1) ) + ...
            -Cuy*( (Iv_uy*cv1+yIv_uy_1).*(Au_uy*u1+yAu_uy_1) );
    kv1_c = -Cvx*( (Iu_vx*cu1+yIu_vx_1).*(Av_vx*v1+yAv_vx_1) ) + ...
            -Cvy*( (Iv_vy*cv1+yIv_vy_1).*(Av_vy*v1+yAv_vy_1) );
    ku2_c = -Cux*( (Iu_ux*cu2+yIu_ux_2).*(Au_ux*u2+yAu_ux_2) ) + ...
            -Cuy*( (Iv_uy*cv2+yIv_uy_2).*(Au_uy*u2+yAu_uy_2) );
    kv2_c = -Cvx*( (Iu_vx*cu2+yIu_vx_2).*(Av_vx*v2+yAv_vx_2) ) + ...
            -Cvy*( (Iv_vy*cv2+yIv_vy_2).*(Av_vy*v2+yAv_vy_2) );


    uh    = uhn + dt*Omu_inv.*( b1*( ku1_d + ku1_c + Fx_1 - y_px_1) + ...
                                b2*( ku2_d + ku2_c + Fx_2 - y_px_2));
    vh    = vhn + dt*Omv_inv.*( b1*( kv1_d + kv1_c + Fy_1 - y_py_1) + ...
                                b2*( kv2_d + kv2_c + Fy_2 - y_py_2));
   
elseif (radau_type==0)
    % stiffly accurate, we don't need additional Poisson solve
    uh    = u2;
    vh    = v2;

end

if (radau_type==1)
    psi1    = a_inv_11*c1*p1 + a_inv_12*c2*p2;
    psi2    = a_inv_21*c1*p1 + a_inv_22*c2*p2;    
    uh = uh - dt*Omu_inv.*( Gx*(b1*psi1+b2*psi2) );
    vh = vh - dt*Omv_inv.*( Gy*(b1*psi1+b2*psi2) );      
end

V    = [uh;vh];

% with the following formulation we don't have to recompute the right hand side with
% the new velocity vector, but simply use the inverse of the Butcher
% tableau; this formulation is however more sensitive to round-off error
% propagation!
% uh = uhn + b1*(a_inv_11*(u1-uhn) + a_inv_12*(u2-uhn)) + b2*(a_inv_21*(u1-uhn) + a_inv_22*(u2-uhn));
% vh = vhn + b1*(a_inv_11*(v1-vhn) + a_inv_12*(v2-vhn)) + b2*(a_inv_21*(v1-vhn) + a_inv_22*(v2-vhn));
% V  = [uh;vh];


% radau_type=0: stiffly accurate, not necessary
% radau_type=1: only necessary for unsteady BC
% radau_type=2: always necessary

if ( (radau_type==1 && BC_unsteady==1) || radau_type==2 ) 
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


end

% this is the standard approach from Hairer '89
% p       = pn + b1* (a_inv_11*(psi1-pn) + a_inv_12*(psi2-pn)) + ...
%                b2* (a_inv_21*(psi1-pn) + a_inv_22*(psi2-pn));

% this is the formula resulting from the approach integral averages ->
% point values; it equals the above formulation for Radau IIA since C(2) holds and 
% furthermore c2=1
p = p1*(2-c2)/(c1-c2) + p2*(c1-2)/(c1-c2);


% store new velocity field in V_ep
V_ep(:,1:method_startup_no-1) = V_ep(:,2:method_startup_no);
V_ep(:,method_startup_no)     = V;
           
           
pressure_additional_solve;