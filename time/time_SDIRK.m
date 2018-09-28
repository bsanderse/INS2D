% SDIRK method of order 3
% solution of saddlepoint system,
% nonlinear_Newton iteration

% solve for U_i

% RK coefficients
% due to a12=0 the stages U_1 and U_2 can be solved sequentially, but this
% has not yet been implemented
gamma  = 1/2 + sqrt(3)/6;
a11    = gamma;
a12    = 0;
a21    = 1-2*gamma;
a22    = gamma;

b1     = 1/2;
b2     = 1/2;

c1     = gamma;
c2     = 1-gamma;

a_inv_11 = 1/gamma;
a_inv_12 = 0;
a_inv_21 = (-1+2*gamma)/(gamma^2);
a_inv_22 = 1/gamma;


f1     = zeros(Nu+Nv,1);
f2     = zeros(Nu+Nv,1);
Z1     = spalloc(Np,Nu+Nv,0);
Z2     = spalloc(2*Np,2*Np,0);

% velocities at time level n
uhn    = uh;
vhn    = vh;
Vn     = V;
pn     = p;
tn     = t;

% initialization
u1     = uh;
v1     = vh;
u2     = uh;
v2     = vh;
p1     = p;
p2     = p;


if (nonlinear_build_matrix==0)
      diag1 = [Omu;Omv;Omu;Omv];
      CD    = spdiags(diag1,0,2*(Nu+Nv),2*(Nu+Nv));
end

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


while (error_nonlinear > nonlinear_acc && i<nonlinear_maxit)
 
      % iterate until u^(n+1) = u^i = u^(i-1); this is different from
      % u^(n)!
    
           
      %% rhs 1

      
      % convective terms, u-component
      uIux       = Iu_ux*u1 + yIu_ux_1; % convective velocity, u_bar
      uAux       = Au_ux*u1 + yAu_ux_1;
      yConv_ux_1 = Cux*( uIux.*uAux );

      vIuy       = Iv_uy*v1 + yIv_uy_1; % convective velocity, v_bar
      uAuy       = Au_uy*u1 + yAu_uy_1;
      yConv_uy_1 = Cuy*( vIuy.*uAuy );

      % convective terms, v-component
      uIvx       = Iu_vx*u1 + yIu_vx_1; % convective velocity, u_bar  
      vAvx       = Av_vx*v1 + yAv_vx_1;
      yConv_vx_1 = Cvx*( uIvx.*vAvx );

      vIvy       = Iv_vy*v1 + yIv_vy_1; % convective velocity, v_bar
      vAvy       = Av_vy*v1 + yAv_vy_1;   
      yConv_vy_1 = Cvy*( vIvy.*vAvy );



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
                    
      
      %% rhs 2
      
      % convective terms, u-component
      uIux       = Iu_ux*u2 + yIu_ux_2; % convective velocity, u_bar
      uAux       = Au_ux*u2 + yAu_ux_2;
      yConv_ux_2 = Cux*( uIux.*uAux );

      vIuy       = Iv_uy*v2 + yIv_uy_2; % convective velocity, v_bar
      uAuy       = Au_uy*u2 + yAu_uy_2;
      yConv_uy_2 = Cuy*( vIuy.*uAuy );

      % convective terms, v-component
      uIvx       = Iu_vx*u2 + yIu_vx_2; % convective velocity, u_bar  
      vAvx       = Av_vx*v2 + yAv_vx_2;
      yConv_vx_2 = Cvx*( uIvx.*vAvx );

      vIvy       = Iv_vy*v2 + yIv_vy_2; % convective velocity, v_bar
      vAvy       = Av_vy*v2 + yAv_vy_2;   
      yConv_vy_2 = Cvy*( vIvy.*vAvy );     

           
      if (nonlinear_build_matrix==1)
          
          C1         = Cux*spdiags(uIux,0,N1,N1);   
          C2         = Cux*spdiags(uAux,0,N1,N1)*nonlinear_Newton;
          Conv_ux_33 = C1*Au_ux + C2*Iu_ux;
          
          C1         = Cuy*spdiags(vIuy,0,N2,N2);   
          C2         = Cuy*spdiags(uAuy,0,N2,N2)*nonlinear_Newton;
          Conv_uy_33 = C1*Au_uy;
          Conv_uy_34 = C2*Iv_uy;   
          
          C1         = Cvx*spdiags(uIvx,0,N3,N3);
          C2         = Cvx*spdiags(vAvx,0,N3,N3)*nonlinear_Newton;
          Conv_vx_43 = C2*Iu_vx;
          Conv_vx_44 = C1*Av_vx;
          
          C1         = Cvy*spdiags(vIvy,0,N4,N4);
          C2         = Cvy*spdiags(vAvy,0,N4,N4)*nonlinear_Newton;       
          Conv_vy_44 = C1*Av_vy + C2*Iv_vy;   
                  
          % assemble matrix 2   
          CD_33      = - Diffu + Conv_ux_33 + Conv_uy_33;
          CD_34      = spalloc(Nu,Nv,0) + Conv_uy_34;
          CD_43      = spalloc(Nv,Nu,0) + Conv_vx_43;
          CD_44      = - Diffv + Conv_vx_44 + Conv_vy_44; 

          CD_2       = [CD_33 CD_34; CD_43 CD_44];
          
      end

      
      % du/dt = k(u), u_i = u^n + dt*(a11*ku1 + a12*ku2)
      ku1        = -yConv_ux_1 - yConv_uy_1 + Diffu*u1 + yDiffu_1 + Fx_1;% - Gx*p1 - y_px_1;
      kv1        = -yConv_vx_1 - yConv_vy_1 + Diffv*v1 + yDiffv_1 + Fy_1;% - Gy*p1 - y_py_1;
      ku2        = -yConv_ux_2 - yConv_uy_2 + Diffu*u2 + yDiffu_2 + Fx_2;% - Gx*p2 - y_px_2;
      kv2        = -yConv_vx_2 - yConv_vy_2 + Diffv*v2 + yDiffv_2 + Fy_2;% - Gy*p2 - y_py_2;
      
      
      % assemble rhs 1
      f1(1:Nu)        =  Omu.*(uhn-u1) + dt*(a11*ku1 + a12*ku2 - c1*(Gx*p1 + y_px_1));                       
      f1(Nu+1:Nu+Nv)  =  Omv.*(vhn-v1) + dt*(a11*kv1 + a12*kv2 - c1*(Gy*p1 + y_py_1));
    
      % assemble rhs 2 
      f2(1:Nu)        =  Omu.*(uhn-u2) + dt*(a21*ku1 + a22*ku2 - c2*(Gx*p2 + y_px_2));                        
      f2(Nu+1:Nu+Nv)  =  Omv.*(vhn-v2) + dt*(a21*kv1 + a22*kv2 - c2*(Gy*p2 + y_py_2));

      %% assemble total right-hand side, this is -1*residual
      fM1       = -M*[u1;v1] - yM1;
      fM2       = -M*[u2;v2] - yM2;

      
      %% total matrix       
      if (nonlinear_build_matrix==1)
          diag1 = [Omu;Omv;Omu;Omv];
          CD    = spdiags(diag1,0,2*(Nu+Nv),2*(Nu+Nv));
          CD    = CD + dt*[a11*CD_1 a12*CD_2; a21*CD_1 a22*CD_2]; 
      end
            
      grad    = [c1*dt*G Z1'; Z1' c2*dt*G];     
      div     = [M Z1; Z1 M];        
      Z       = [CD grad; div Z2];        
      f       = [f1;f2;fM1;fM2];

      % solve
      dq        = Z\f;

      du1       = dq(1:Nu);
      dv1       = dq(Nu+1:Nu+Nv);
      du2       = dq(Nu+Nv+1:2*Nu+Nv);
      dv2       = dq(2*Nu+Nv+1:2*(Nu+Nv));
      dp1       = dq(2*(Nu+Nv)+1:2*(Nu+Nv)+Np);
      dp2       = dq(2*(Nu+Nv)+Np+1:2*(Nu+Nv)+2*Np);


      u1        = u1 + du1;
      v1        = v1 + dv1;
      u2        = u2 + du2;
      v2        = v2 + dv2; 
      p1        = p1 + dp1;
      p2        = p2 + dp2;        
        
 
      % note that dp does not necessarily go to zero
      % divide the residual by a norm based on rhs ?
      k_norm = [Om.*Vn;Om.*Vn];
%       error_nonlinear = norm(f,inf) /norm(k_norm,inf);
      error_nonlinear = norm(f,inf);
%       error_nonlinear = max(norm(f),norm(dq))/norm(k_norm)
      i = i+1;
      
end
% pause;
if (i>= nonlinear_maxit || error_nonlinear>nonlinear_acc || isnan(error_nonlinear))
    error(['nonlinear iterations did not converge in ' num2str(nonlinear_maxit) ' iterations']);
end

% update with F1 and F2 
% uh   = uhn + dt*Omu_inv.*(b1*ku1 + b2*ku2);
% vh   = vhn + dt*Omv_inv.*(b1*kv1 + b2*kv2);
% V    = [uh;vh];

% with the following formulation we don't have to recompute the right hand side with
% the new velocity vector, but simply use the inverse of the Butcher
% tableau
uh = uhn + b1*(a_inv_11*(u1-uhn) + a_inv_12*(u2-uhn)) + b2*(a_inv_21*(u1-uhn) + a_inv_22*(u2-uhn));
vh = vhn + b1*(a_inv_11*(v1-vhn) + a_inv_12*(v2-vhn)) + b2*(a_inv_21*(v1-vhn) + a_inv_22*(v2-vhn));
V  = [uh;vh];


p1_temp = a_inv_11*(p1 - pn) + a_inv_12*(p2 - pn);
p2_temp = a_inv_21*(p1 - pn) + a_inv_22*(p2 - pn);
p       = pn + b1*p1_temp + b2*p2_temp;

% if (BC_unsteady == 1)
%     % make V satisfy the incompressibility constraint at n+1
%     t = tn + dt;
%     boundary_conditions;
%     interpolate_bc;
%     operator_bc_divergence;
%     t = tn;
% 
%     f       = (1/dt)*(M*V + yM);
% 
%     pressure_poisson;
% 
%     V       = V - dt*Om_inv.*(G*dp);
%     uh      = V(1:Nu);
%     vh      = V(Nu+1:Nu+Nv);
% 
%     p       = pn + dp;
% 
% end

% store velocity fields if this method is used as startup 
% for linear extrapolation
% if (method_temp == 92)
%     V_ep(:,n) = V;  
% end

pressure_additional_solve;