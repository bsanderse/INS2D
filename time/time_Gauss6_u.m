% fully (energy) conservative time stepping
% solution of saddlepoint system,
% nonlinear_Newton iteration

% solve for U_i

% reversible in time (inviscid flow)

% RK coefficients
a11    = 5/36;
a12    = 2/9 - sqrt(15)/15;
a13    = 5/36 - sqrt(15)/30;
a21    = 5/36 + sqrt(15)/24;
a22    = 2/9;
a23    = 5/36 - sqrt(15)/24;
a31    = 5/36 + sqrt(15)/30;
a32    = 2/9 + sqrt(15)/15;
a33    = 5/36;


b1     = 5/18;
b2     = 4/9;
b3     = 5/18;

c1     = 1/2 - sqrt(15)/10;
c2     = 1/2;
c3     = 1/2 + sqrt(15)/10;

a_inv_11 = 5; 
a_inv_12 = (4/3)*sqrt(15)-4;
a_inv_13 = -(4/3)*sqrt(15)+5;
a_inv_21 = -(5/6)*sqrt(15)-5/2;
a_inv_22 = 2;
a_inv_23 = (5/6)*sqrt(15)-5/2;
a_inv_31 = (4/3)*sqrt(15)+5;
a_inv_32 = -(4/3)*sqrt(15)-4;
a_inv_33 = 5;


f1     = zeros(Nu+Nv,1);
f2     = zeros(Nu+Nv,1);
f3     = zeros(Nu+Nv,1);
Z1     = spalloc(Np,Nu+Nv,0);
Z2     = spalloc(3*Np,3*Np,0);

if (n==2 || nonlinear_startingvalues==0)
    % initialization
    u1     = uh;
    v1     = vh;
    u2     = uh;
    v2     = vh;
    u3     = uh;
    v3     = vh;
    p1     = p;
    p2     = p;
    p3     = p;
else
    % better initialization using values from previous stage
    s        = 3;
    c_list   = [0 c1 c2 c3 1];
    id       = eye(s+2);
    for is=1:s+2
        poly_coeff     = polyfit(c_list,id(is,:),s+1);
        coeff(is,1:s)  = polyval(poly_coeff,1+[c1 c2 c3]);
    end
    u1_temp     = coeff(1,1)*uhn + coeff(2,1)*u1 + coeff(3,1)*u2 + coeff(4,1)*u3 + coeff(5,1)*uh;
    u2_temp     = coeff(1,2)*uhn + coeff(2,2)*u1 + coeff(3,2)*u2 + coeff(4,2)*u3 + coeff(5,2)*uh;
    u3_temp     = coeff(1,3)*uhn + coeff(2,3)*u1 + coeff(3,3)*u2 + coeff(4,3)*u3 + coeff(5,3)*uh;
    v1_temp     = coeff(1,1)*vhn + coeff(2,1)*v1 + coeff(3,1)*v2 + coeff(4,1)*v3 + coeff(5,1)*vh;
    v2_temp     = coeff(1,2)*vhn + coeff(2,2)*v1 + coeff(3,2)*v2 + coeff(4,2)*v3 + coeff(5,2)*vh;
    v3_temp     = coeff(1,3)*vhn + coeff(2,3)*v1 + coeff(3,3)*v2 + coeff(4,3)*v3 + coeff(5,3)*vh;
    p1_temp     = coeff(1,1)*pn  + coeff(2,1)*p1 + coeff(3,1)*p2 + coeff(4,1)*p3 + coeff(5,1)*p;
    p2_temp     = coeff(1,2)*pn  + coeff(2,2)*p1 + coeff(3,2)*p2 + coeff(4,2)*p3 + coeff(5,2)*p;    
    p3_temp     = coeff(1,3)*pn  + coeff(2,3)*p1 + coeff(3,3)*p2 + coeff(4,3)*p3 + coeff(5,3)*p;    
    u1 = u1_temp;
    u2 = u2_temp;
    u3 = u3_temp;
    v1 = v1_temp;
    v2 = v2_temp;
    v3 = v3_temp;
    p1 = p1_temp;
    p2 = p2_temp;
    p3 = p3_temp;
end

if (nonlinear_build_matrix==0)
      diag1 = [Omu;Omv;Omu;Omv;Omu;Omv];
      CD    = spdiags(diag1,0,3*(Nu+Nv),3*(Nu+Nv));
end

% velocities at time level n
uhn    = uh;
vhn    = vh;
Vn     = V;
pn     = p;
tn     = t;

error_nonlinear    = zeros(nonlinear_maxit,1);
error_nonlinear(1) = 1;% iterate to remove linearization error; necessary to obtain full
% conservation and reversibility; if not removed still fully conservative
% but not reversible, and maybe not 2nd order
i=0;


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


t = tn + c3*dt;
boundary_conditions;
interpolate_bc;
operator_bc_divergence;
operator_bc_momentum;
force;
t = tn;

yIu_ux_3   = yIu_ux;
yAu_ux_3   = yAu_ux;
yIv_uy_3   = yIv_uy;
yAu_uy_3   = yAu_uy;
yIu_vx_3   = yIu_vx;
yAv_vx_3   = yAv_vx;
yIv_vy_3   = yIv_vy;
yAv_vy_3   = yAv_vy;
yDiffu_3   = yDiffu;
yDiffv_3   = yDiffv;
Fx_3       = Fx;
Fy_3       = Fy;
y_px_3     = y_px;
y_py_3     = y_py;
yM3        = yM;   


while (i<nonlinear_maxit)
 
      i = i+1;
      
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
                  
          % assemble matrix 2   
          CD_11      = - Diffu + Conv_ux_11 + Conv_uy_11;
          CD_12      = spalloc(Nu,Nv,0) + Conv_uy_12;
          CD_21      = spalloc(Nv,Nu,0) + Conv_vx_21;
          CD_22      = - Diffv + Conv_vx_22 + Conv_vy_22; 

          CD_2       = [CD_11 CD_12; CD_21 CD_22];
          
      end

      %% rhs 3
      
      % convective terms, u-component
      uIux       = Iu_ux*u3 + yIu_ux_3; % convective velocity, u_bar
      uAux       = Au_ux*u3 + yAu_ux_3;
      yConv_ux_3 = Cux*( uIux.*uAux );

      vIuy       = Iv_uy*v3 + yIv_uy_3; % convective velocity, v_bar
      uAuy       = Au_uy*u3 + yAu_uy_3;
      yConv_uy_3 = Cuy*( vIuy.*uAuy );

      % convective terms, v-component
      uIvx       = Iu_vx*u3 + yIu_vx_3; % convective velocity, u_bar  
      vAvx       = Av_vx*v3 + yAv_vx_3;
      yConv_vx_3 = Cvx*( uIvx.*vAvx );

      vIvy       = Iv_vy*v3 + yIv_vy_3; % convective velocity, v_bar
      vAvy       = Av_vy*v3 + yAv_vy_3;   
      yConv_vy_3 = Cvy*( vIvy.*vAvy );     

           
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
                  
          % assemble matrix 2   
          CD_11      = - Diffu + Conv_ux_11 + Conv_uy_11;
          CD_12      = spalloc(Nu,Nv,0) + Conv_uy_12;
          CD_21      = spalloc(Nv,Nu,0) + Conv_vx_21;
          CD_22      = - Diffv + Conv_vx_22 + Conv_vy_22; 

          CD_3       = [CD_11 CD_12; CD_21 CD_22];
          
      end
      
      % du/dt = k(u), u_i = u^n + dt*(a11*ku1 + a12*ku2)
      ku1        = -yConv_ux_1 - yConv_uy_1 + Diffu*u1 + yDiffu_1 + Fx_1;% - Gx*p1 - y_px_1;
      kv1        = -yConv_vx_1 - yConv_vy_1 + Diffv*v1 + yDiffv_1 + Fy_1;% - Gy*p1 - y_py_1;
      ku2        = -yConv_ux_2 - yConv_uy_2 + Diffu*u2 + yDiffu_2 + Fx_2;% - Gx*p2 - y_px_2;
      kv2        = -yConv_vx_2 - yConv_vy_2 + Diffv*v2 + yDiffv_2 + Fy_2;% - Gy*p2 - y_py_2;
      ku3        = -yConv_ux_3 - yConv_uy_3 + Diffu*u3 + yDiffu_3 + Fx_3;% - Gx*p2 - y_px_2;
      kv3        = -yConv_vx_3 - yConv_vy_3 + Diffv*v3 + yDiffv_3 + Fy_3;% - Gy*p2 - y_py_2;      
      
      % assemble rhs 1
      f1(1:Nu)        =  Omu.*(uhn-u1) + dt*(a11*ku1 + a12*ku2 + a13*ku3 - c1*(Gx*p1 + y_px_1));                       
      f1(Nu+1:Nu+Nv)  =  Omv.*(vhn-v1) + dt*(a11*kv1 + a12*kv2 + a13*kv3 - c1*(Gy*p1 + y_py_1));
    
      % assemble rhs 2 
      f2(1:Nu)        =  Omu.*(uhn-u2) + dt*(a21*ku1 + a22*ku2 + a23*ku3 - c2*(Gx*p2 + y_px_2));                        
      f2(Nu+1:Nu+Nv)  =  Omv.*(vhn-v2) + dt*(a21*kv1 + a22*kv2 + a23*kv3 - c2*(Gy*p2 + y_py_2));

      % assemble rhs 2 
      f3(1:Nu)        =  Omu.*(uhn-u3) + dt*(a31*ku1 + a32*ku2 + a33*ku3 - c3*(Gx*p3 + y_px_3));                        
      f3(Nu+1:Nu+Nv)  =  Omv.*(vhn-v3) + dt*(a31*kv1 + a32*kv2 + a33*kv3 - c3*(Gy*p3 + y_py_3));      
      
      
      %% assemble total right-hand side, this is -1*residual
      fM1       = -M*[u1;v1] - yM1;
      fM2       = -M*[u2;v2] - yM2;
      fM3       = -M*[u3;v3] - yM3;

      f       = [f1;f2;f3;fM1;fM2;fM3];
      
      error_nonlinear(i) = norm(f,inf);
      if (error_nonlinear(i) <= nonlinear_acc || ...
          error_nonlinear(i)/error_nonlinear(1)<= nonlinear_relacc)
          break;
      end   
      
      %% total matrix       
      if (nonlinear_build_matrix==1)
          diag1 = [Omu;Omv;Omu;Omv;Omu;Omv];
          CD    = spdiags(diag1,0,3*(Nu+Nv),3*(Nu+Nv));
          CD    = CD + dt*[a11*CD_1 a12*CD_2 a13*CD_3; ...
                           a21*CD_1 a22*CD_2 a23*CD_3; ...
                           a31*CD_1 a32*CD_2 a33*CD_3]; 
      end
            
      grad    = [c1*dt*G Z1' Z1'; Z1' c2*dt*G Z1'; Z1' Z1' c3*dt*G];     
      div     = [M Z1 Z1; Z1 M Z1; Z1 Z1 M]; 
      Z       = [CD grad; div Z2];        

      % solve
      dq        = Z\f;

      du1       = dq(1:Nu);
      dv1       = dq(Nu+1:Nu+Nv);
      du2       = dq(Nu+Nv+1:2*Nu+Nv);
      dv2       = dq(2*Nu+Nv+1:2*(Nu+Nv));
      du3       = dq(2*(Nu+Nv)+1:3*Nu+2*Nv);
      dv3       = dq(3*Nu+2*Nv+1:3*(Nu+Nv));      
      dp1       = dq(3*(Nu+Nv)+1:3*(Nu+Nv)+Np);
      dp2       = dq(3*(Nu+Nv)+Np+1:3*(Nu+Nv)+2*Np);
      dp3       = dq(3*(Nu+Nv)+2*Np+1:3*(Nu+Nv)+3*Np);


      u1        = u1 + du1;
      v1        = v1 + dv1;
      u2        = u2 + du2;
      v2        = v2 + dv2; 
      u3        = u3 + du3;
      v3        = v3 + dv3;       
      p1        = p1 + dp1;
      p2        = p2 + dp2;   
      p3        = p3 + dp3; 
        
 
      % note that dp does not necessarily go to zero
      % divide the residual by a norm based on rhs ?
%       k_norm = [Om.*Vn;Om.*Vn;Om.*Vn];
%       error_nonlinear = norm(f,inf) /norm(k_norm,inf);
%       error_nonlinear = norm(f,inf);
%       error_nonlinear = max(norm(f),norm(dq))/norm(k_norm)
%       i = i+1;
      
end
if ((error_nonlinear(i)>nonlinear_acc && error_nonlinear(i)/error_nonlinear(1)> nonlinear_relacc) ...
     || i>= nonlinear_maxit || isnan(error_nonlinear(i)))
    error(['nonlinear iterations did not converge in ' num2str(nonlinear_maxit) ' iterations']);
    error_nonlinear(i)
    error_nonlinear(i)/error_nonlinear(1)
end
nonlinear_its(n) = i-1;

% update with F1 and F2 
% uh   = uhn + dt*Omu_inv.*(b1*ku1 + b2*ku2);
% vh   = vhn + dt*Omv_inv.*(b1*kv1 + b2*kv2);
% V    = [uh;vh];

% with the following formulation we don't have to recompute the right hand side with
% the new velocity vector, but simply use the inverse of the Butcher
% tableau
uh = uhn + b1*(a_inv_11*(u1-uhn) + a_inv_12*(u2-uhn) + a_inv_13*(u3-uhn)) + ...
           b2*(a_inv_21*(u1-uhn) + a_inv_22*(u2-uhn) + a_inv_23*(u3-uhn)) + ...
           b3*(a_inv_31*(u1-uhn) + a_inv_32*(u2-uhn) + a_inv_33*(u3-uhn));
vh = vhn + b1*(a_inv_11*(v1-vhn) + a_inv_12*(v2-vhn) + a_inv_13*(v3-vhn)) + ...
           b2*(a_inv_21*(v1-vhn) + a_inv_22*(v2-vhn) + a_inv_23*(v3-vhn)) + ...
           b3*(a_inv_31*(v1-vhn) + a_inv_32*(v2-vhn) + a_inv_33*(v3-vhn));
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


% p1_temp = a_inv_11*c1*p1 + a_inv_12*c2*p2 + a_inv_13*c3*p3;
% p2_temp = a_inv_21*c1*p1 + a_inv_22*c2*p2 + a_inv_23*c3*p3;
% p3_temp = a_inv_31*c1*p1 + a_inv_32*c2*p2 + a_inv_33*c3*p3;

% using a Lagrange polynomial through p1, p2 and p3 (Hairer p. 505)
% leads to an s-th order accurate pressure
% p = p1_temp*(1-c2)*(1-c3)/((c1-c2)*(c1-c3)) + ...
%     p2_temp*(1-c1)*(1-c3)/((c2-c1)*(c2-c3)) + ...
%     p3_temp*(1-c1)*(1-c2)/((c3-c1)*(c3-c2));

% this is equivalent to our approach,
p = p1*(3-2*c3-2*c2+c2*c3)/((c1-c2)*(c1-c3)) + ...
   -p2*(3-2*c3-2*c1+c1*c3)/((c1-c2)*(c2-c3)) + ...
    p3*(3-2*c2-2*c1+c1*c2)/((c1-c3)*(c2-c3));


pressure_additional_solve;