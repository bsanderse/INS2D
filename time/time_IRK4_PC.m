% fully (energy) conservative time stepping
% solution of saddlepoint system,
% nonlinear_Newton iteration

% reversible in time (inviscid flow)

% RK coefficients
a11    = 1/4;
a12    = 1/4-sqrt(3)/6;
a21    = 1/4+sqrt(3)/6;
a22    = 1/4;

b1     = 1/2;
b2     = 1/2;

c1     = 1/2-sqrt(3)/6;
c2     = 1/2+sqrt(3)/6;

a_inv_11 = 3;
a_inv_12 = -3+2*sqrt(3);
a_inv_21 = -3-2*sqrt(3);
a_inv_22 = 3;

% CN     = 0.5;
f      = zeros(2*(Nu+Nv+Np),1);
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
ku1    = u1;
kv1    = v1;
ku2    = u2;
kv2    = v2;
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
operator_boundary_conditions;
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
operator_boundary_conditions;
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


while (error_nonlinear > nonlinear_acc  && i<nonlinear_maxit)
 
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

      % assemble rhs 1
      f1(1:Nu)   =  Omu.*ku1 + yConv_ux_1 + yConv_uy_1 + ...
                   -Diffu*u1 - yDiffu_1 - Fx_1 + Gx*p1 + y_px_1;
                         
      f1(Nu+1:Nu+Nv)  =  Omv.*kv1 + yConv_vx_1 + yConv_vy_1 + ...
                        -Diffv*v1 - yDiffv_1 - Fy_1 + Gy*p1 + y_py_1;

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

          % BLOCKS 13 - 14 - 23 - 24
          % same blocks as 11-12-21-22

          % assemble matrix 1
          
          CD_11      = - Diffu + Conv_ux_11 + Conv_uy_11;
          CD_12      = spalloc(Nu,Nv,0) + Conv_uy_12;
          CD_21      = spalloc(Nv,Nu,0) + Conv_vx_21;
          CD_22      = - Diffv  + Conv_vx_22 + Conv_vy_22; 

          CD_1       = [CD_11 CD_12; CD_21 CD_22];
% 
%           CD_13      = a12*dt*(- Diffu + Conv_ux_11 + Conv_uy_11); 
%           CD_14      = a12*dt*Conv_uy_12;
%           CD_23      = a12*dt*Conv_vx_21;
%           CD_24      = a12*dt*(- Diffv + Conv_vx_22 + Conv_vy_22); 

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

    
      % assemble rhs 2 
      f2(1:Nu)   =  Omu.*ku2 + yConv_ux_2 + yConv_uy_2 + ...
                   -Diffu*u2 - yDiffu_2 - Fx_2 + Gx*p2 + y_px_2;
               
      f2(Nu+1:Nu+Nv)  =  Omv.*kv2 + yConv_vx_2 + yConv_vy_2 + ...
                        -Diffv*v2 - yDiffv_2 - Fy_2 + Gy*p2 + y_py_2;
           
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
          
          % BLOCKS 31 - 32 - 41 - 42
          % same blocks as 33 - 34 - 43 - 44   
          
          % assemble matrix 2   
          CD_33      = - Diffu + Conv_ux_33 + Conv_uy_33;
          CD_34      = spalloc(Nu,Nv,0) + Conv_uy_34;
          CD_43      = spalloc(Nv,Nu,0) + Conv_vx_43;
          CD_44      = - Diffv + Conv_vx_44 + Conv_vy_44; 

          CD_2       = [CD_33 CD_34; CD_43 CD_44];
          
%           CD_31      = a21*dt*(- Diffu + Conv_ux_33 + Conv_uy_33);
%           CD_32      = a21*dt*Conv_uy_34;
%           CD_41      = a21*dt*Conv_vx_43;
%           CD_42      = a21*dt*(- Diffv + Conv_vx_44 + Conv_vy_44);            
          
      end
               

      %% total matrix       
      if (nonlinear_build_matrix==1)
          diag1 = [Omu;Omv;Omu;Omv];
          CD    = spdiags(diag1,0,2*(Nu+Nv),2*(Nu+Nv));
          CD    = CD + dt*[a11*CD_1 a12*CD_1; a21*CD_2 a22*CD_2]; 
      end
      
%       grad      = [G Z1'; Z1' G];      
%       div       = [M Z1; Z1 M];
%       Z         = [CD grad; div Z2]; 
%       if (i==1)
%           D = eigs(Z)
%       end

      %% assemble total right-hand side, this is -1*residual
%       fM1       = M*[ku1;kv1] + yM1;
%       fM2       = M*[ku2;kv2] + yM2;
      f         = [f1;f2];%fM1;fM2];

      % solve
      dq        = -CD\f;
 
      dku1      = dq(1:Nu);
      dkv1      = dq(Nu+1:Nu+Nv);
      dku2      = dq(Nu+Nv+1:2*Nu+Nv);
      dkv2      = dq(2*Nu+Nv+1:2*(Nu+Nv));
%       dp1       = dq(2*(Nu+Nv)+1:2*(Nu+Nv)+Np);
%       dp2       = dq(2*(Nu+Nv)+Np+1:2*(Nu+Nv)+2*Np);
      

      ku1       = ku1 + dku1;
      kv1       = kv1 + dkv1;
      ku2       = ku2 + dku2;
      kv2       = kv2 + dkv2; 
%       p1        = p1 + dp1;
%       p2        = p2 + dp2;
      
      u1        = uhn + dt*(a11*ku1 + a12*ku2);
      v1        = vhn + dt*(a11*kv1 + a12*kv2);
      u2        = uhn + dt*(a21*ku1 + a22*ku2);
      v2        = vhn + dt*(a21*kv1 + a22*kv2);       
      
      % note that dp does not necessarily go to zero
      % divide the residual by a norm based on rhs ?
      k_norm = [Om.*Vn;Om.*Vn];
%       error_nonlinear = norm(f,inf) /norm(k_norm,inf);
      error_nonlinear = norm(f,inf)
%       error_nonlinear = max(norm(f),norm(dq))/norm(k_norm)
      i = i+1;
%       pause;
      
end

if (i>= nonlinear_maxit || error_nonlinear>nonlinear_acc || isnan(error_nonlinear))
    error(['nonlinear iterations did not converge in ' num2str(nonlinear_maxit) ' iterations']);
end

uh   = uhn + dt*(b1*ku1 + b2*ku2);
vh   = vhn + dt*(b1*kv1 + b2*kv2);
R    = [uh;vh];

pressure_correction;
p_old = p;
p     = p + dp;

uh_old = uh;
vh_old = vh;
uh     = V(1:Nu);
vh     = V(Nu+1:Nu+Nv);
% 
% p1_temp = a_inv_11*(p1 - pn) + a_inv_12*(p2 - pn);
% p2_temp = a_inv_21*(p1 - pn) + a_inv_22*(p2 - pn);
% p       = pn + b1*p1_temp + b2*p2_temp;

%%
% make the pressure at n+1 compatible with the velocity field. this should
% also result in 4th order pressure

if (p_add_solve==1)
    % evaluate BC and force at final time    
    t = tn + dt;
    boundary_conditions;
    operator_boundary_conditions;
    force;
    t = tn;    
    
    % additional poisson solve for 4th order accurate pressure
    % convection
    cu = uh;
    cv = vh;
    convection;

    % diffusion
    d2u    = Diffu*uh + yDiffu;
    d2v    = Diffv*vh + yDiffv;


    Ru =  - du2dx - duvdy + d2u + Fx - y_px;    
    Rv =  - duvdx - dv2dy + d2v + Fy - y_py;
    R  = Om_inv.*[Ru;Rv];
    f  = (M*R + yM);

    dp = zeros(Np,1);

    if (poisson==1)
        % using pre-determined LU decomposition
        b   = L\f;
        dp  = U\b;
    elseif (poisson==2)
      % using preconditioned CG
        t_pressure = toc;
        [dp, flag, norm1, iter]  = pcg(-A,-f,CG_acc,CG_maxit,L,L',dp);
        t_solve    = toc - t_pressure;
    elseif (poisson==3)  
        t_pressure = toc;
        [dp,iter,norm1,norm2]=cg(B,int64(dia),int64(ndia),f,CG_acc,int64(Np),dp,int64(CG_maxit));
        t_solve    = toc - t_pressure;    
    elseif (poisson==4)
        t_pressure = toc;  
        [dp,iter,norm1,norm2]=cg_matlab(A,f,CG_acc,Np,dp,CG_maxit,A_pc);
        t_solve    = toc - t_pressure;
    elseif (poisson==5)
        t_pressure = toc;
        PetscBinaryWrite(PS,-f');
        t_write    = toc - t_pressure;

        t_pressure = toc;   
        dp     = PetscBinaryRead(PS);
        iter   = PetscBinaryRead(PS);
        norm1  = PetscBinaryRead(PS);   
        t_read = toc - t_pressure;

        dp     = dp';
        iter   = iter(1);
        norm1  = norm1(1);  
    end

    % write information on pressure solve
    if (poisson==2 || poisson==3 || poisson==4)
        fprintf(fpres,'%-10i %-10i %16.8e %16.8e\n',n,iter,norm1,t_solve);
    elseif (poisson == 5)
        fprintf(fpres,'%-10i %-10i %16.8e %16.8e %16.8e\n',n,iter,norm1,t_write,t_read);
    end    

    p = dp;
end