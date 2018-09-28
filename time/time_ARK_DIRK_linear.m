% Radau IIA method of order 3
% solution of saddlepoint system,
% nonlinear_Newton iteration

% solve for U_i

dirk_type = type_list(jj); % 0: diffusion, 1: convection, 2: mix


% RK coefficients for convection
if (dirk_type == 0)
    a11c    = 1/4;
    a12c    = 0;
    a21c    = 5/12;
    a22c    = 1/3;
else
    a11c    = 1/4;
    a12c    = 0;
    a21c    = 1/2;
    a22c    = 1/4;
end
    
% RK coefficients for diffusion
if (dirk_type == 1)
    a11d    = 1/4;
    a12d    = 0;
    a21d    = 1/2;
    a22d    = 1/4;
else
    a11d    = 1/4;
    a12d    = 0;
    a21d    = 5/12;
    a22d    = 1/3;
end

% c-coefficients are same for convection and diffusion
c1     = 1/4;
c2     = 3/4;

b1     = 1/2;
b2     = 1/2;

% inverse of A-coefficients, used for pressure
if (dirk_type == 1)
    % convection:
    a_inv_11 = 4;
    a_inv_12 = 0;
    a_inv_21 = -8;
    a_inv_22 = 4;
else
    % diffusion:
    a_inv_11 = 4;
    a_inv_12 = 0;
    a_inv_21 = -5;
    a_inv_22 = 3;
end

% coefficients for extrapolation of convective velocity
% Vandermonde matrix, in case of constant time step
% we could precompute the coefficients at the start of the time loop
VM_i    = method_startup_no-1:-1:0; 
VM_i    = VM_i/2; % only for DIRK methods: we use the intermediate velocity field Vi in the update

VM1_row = c1+VM_i;
% VM2_row = c2+VM_i;
VM1     = vander(VM1_row)';
% VM2     = vander(VM2_row)';
VM_rhs  = zeros(method_startup_no,1);
VM_rhs(end) = 1;

ep1     = VM1\VM_rhs;
% ep2     = VM2\VM_rhs;

cV1    = V_ep*ep1;
% cV2    = V_ep*ep2;
cu1    = cV1(1:Nu);
cv1    = cV1(Nu+1:Nu+Nv);
% cu2    = cV2(1:Nu);
% cv2    = cV2(Nu+1:Nu+Nv);


% initialization
f1     = zeros(Nu+Nv,1);
f2     = zeros(Nu+Nv,1);
Z1     = spalloc(Np,Nu+Nv,0);
Z2     = spalloc(Np,Np,0);



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

    
  % du/dt = k(u), u_i = u^n + dt*(a11*ku1 + a12*ku2)
  ku1_c        = -yConv_ux_1 - yConv_uy_1;
  kv1_c        = -yConv_vx_1 - yConv_vy_1;
  ku1_d        = yDiffu_1 + Fx_1;
  kv1_d        = yDiffv_1 + Fy_1;


  % assemble rhs 1
  f1(1:Nu)        =  Omu.*uhn + dt*(a11c*ku1_c + a11d*ku1_d ...
                                         - c1*y_px_1);
  f1(Nu+1:Nu+Nv)  =  Omv.*vhn + dt*(a11c*kv1_c + a11d*kv1_d ...
                                         - c1*y_py_1);


   %% assemble total right-hand side, this is -1*residual
   fM1       = - yM1;

      
   %% total matrix       
   if (nonlinear_build_matrix==1)
       diag1 = [Omu;Omv];
       CD    = spdiags(diag1,0,Nu+Nv,Nu+Nv);
       CD    = CD + dt*(a11c*C_1 - a11d*D_1);
   end

   Z       = [CD c1*dt*G; M Z2];        
   f       = [f1;fM1];

   % solve
   dq        = Z\f;

   u1       = dq(1:Nu);
   v1       = dq(Nu+1:Nu+Nv);
   p1       = dq(Nu+Nv+1:Nu+Nv+Np);

   
% F1 based on V1
ku1_d = Diffu*u1 + yDiffu_1;
kv1_d = Diffv*v1 + yDiffv_1;

ku1_c = -Cux*( (Iu_ux*cu1+yIu_ux_1).*(Au_ux*u1+yAu_ux_1) ) + ...
        -Cuy*( (Iv_uy*cv1+yIv_uy_1).*(Au_uy*u1+yAu_uy_1) );
kv1_c = -Cvx*( (Iu_vx*cu1+yIu_vx_1).*(Av_vx*v1+yAv_vx_1) ) + ...
        -Cvy*( (Iv_vy*cv1+yIv_vy_1).*(Av_vy*v1+yAv_vy_1) );

% intermediate second order velocity field:
ui   = uhn + dt*Omu_inv.*( b1*( ku1_d + ku1_c + Fx_1 - y_px_1) - b1*Gx*p1);
vi   = vhn + dt*Omv_inv.*( b1*( kv1_d + kv1_c + Fy_1 - y_py_1) - b1*Gy*p1);
Vi   = [ui; vi];

if (BC_unsteady == 1 || dirk_type==2)
    % make V satisfy the incompressibility constraint at n+1
    t = tn + 2*c1*dt;
    boundary_conditions;
    interpolate_bc;
    operator_bc_divergence;
    t = tn;

    f       = (1/dt)*(M*Vi + yM);

    pressure_poisson;

    Vi       = Vi - dt*Om_inv.*(G*dp);

end

V_ep(:,1) = V_ep(:,2);
V_ep(:,2) = Vi;   

cV2    = V_ep*ep1;
cu2    = cV2(1:Nu);
cv2    = cV2(Nu+1:Nu+Nv);

% next stage
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
  Conv_ux_11 = C1*Au_ux;

  C1         = Cuy*spdiags(vIuy,0,N2,N2);   
  Conv_uy_11 = C1*Au_uy;

  C1         = Cvx*spdiags(uIvx,0,N3,N3);
  Conv_vx_22 = C1*Av_vx;

  C1         = Cvy*spdiags(vIvy,0,N4,N4);
  Conv_vy_22 = C1*Av_vy;  

  
  % assemble matrix

  C_11      = Conv_ux_11 + Conv_uy_11;
  C_12      = spalloc(Nu,Nv,0);
  C_21      = spalloc(Nv,Nu,0);
  C_22      = Conv_vx_22 + Conv_vy_22; 

  C_1       = [C_11 C_12; C_21 C_22];
  D_1       = [Diffu spalloc(Nu,Nv,0); spalloc(Nv,Nu,0) Diffv];
  
      
  % du/dt = k(u), u_i = u^n + dt*(a11*ku1 + a12*ku2)
  ku2_c        = -yConv_ux_2 - yConv_uy_2;
  kv2_c        = -yConv_vx_2 - yConv_vy_2;
  ku2_d        = yDiffu_2 + Fx_2;
  kv2_d        = yDiffv_2 + Fy_2;
      
        
      % assemble rhs 2 
      f2(1:Nu)        =  Omu.*uhn + dt*(a21c*ku1_c + a22c*ku2_c + ...
                                        a21d*ku1_d + a22d*ku2_d + ...
                                             - c2*y_px_2);                        
      f2(Nu+1:Nu+Nv)  =  Omv.*vhn + dt*(a21c*kv1_c + a22c*kv2_c + ...
                                        a21d*kv1_d + a22d*kv2_d + ...
                                             - c2*y_py_2);                                         

      %% assemble total right-hand side, this is -1*residual
      fM2       = - yM2;

      
      %% total matrix       
      if (nonlinear_build_matrix==1)
          diag1 = [Omu;Omv];
          CD    = spdiags(diag1,0,Nu+Nv,Nu+Nv);
          CD    = CD + dt*(a22c*C_1 - a22d*D_1);
      end
            
      Z       = [CD c2*dt*G; M Z2];        
      f       = [f2;fM2];

      % solve
      dq        = Z\f;

      u2       = dq(1:Nu);
      v2       = dq(Nu+1:Nu+Nv);
      p2       = dq(Nu+Nv+1:Nu+Nv+Np);


% F2 based on V2
ku2_d = Diffu*u2 + yDiffu_2;
kv2_d = Diffv*v2 + yDiffv_2;

ku2_c = -Cux*( (Iu_ux*cu2+yIu_ux_2).*(Au_ux*u2+yAu_ux_2) ) + ...
        -Cuy*( (Iv_uy*cv2+yIv_uy_2).*(Au_uy*u2+yAu_uy_2) );
kv2_c = -Cvx*( (Iu_vx*cu2+yIu_vx_2).*(Av_vx*v2+yAv_vx_2) ) + ...
        -Cvy*( (Iv_vy*cv2+yIv_vy_2).*(Av_vy*v2+yAv_vy_2) );

uh    = uhn + dt*Omu_inv.*( b1*( ku1_d + ku1_c + Fx_1 - y_px_1) + ...
                            b2*( ku2_d + ku2_c + Fx_2 - y_px_2));
vh    = vhn + dt*Omv_inv.*( b1*( kv1_d + kv1_c + Fy_1 - y_py_1) + ...
                            b2*( kv2_d + kv2_c + Fy_2 - y_py_2));

% type 0 is NOT stiffly accurate (in contrast to the Radau and Lobatto
% methods); for type 0 and 1 we can avoid an additional projection in case
% of steady BC.
if (dirk_type==0 || dirk_type==1)
    psi1    = a_inv_11*c1*p1 + a_inv_12*c2*p2;
    psi2    = a_inv_21*c1*p1 + a_inv_22*c2*p2;    
    uh = uh - dt*Omu_inv.*( Gx*(b1*psi1+b2*psi2) );
    vh = vh - dt*Omv_inv.*( Gy*(b1*psi1+b2*psi2) );      
end

V    = [uh;vh];


if (BC_unsteady == 1 || dirk_type==2)
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


p = p2;

% store new velocity field in V_ep
V_ep(:,1) = V_ep(:,2);
V_ep(:,2) = V;           


pressure_additional_solve;