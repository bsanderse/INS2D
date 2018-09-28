% fully (energy) conservative time stepping
% solution of saddlepoint system,
% nonlinear_Newton iteration

% not reversible in time (inviscid flow)

% RK coefficients
a11    = 1/2;
b1     = 1;
c1     = 1/2;

a_inv_11 = 2;

f1     = zeros(Nu+Nv,1);
Z2     = spalloc(Np,Np,0);

% velocities at time level n
uhn    = uh;
vhn    = vh;
Vn     = V;
pn     = p;
tn     = t;

% estimate of the convecting velocity at n+1/2
% cu and cv should be the same for time n and time n+1 to conserve energy !!
if (EP==0)
    % Picard
    cu     = uh;
    cv     = vh;
elseif (EP==1)
    % extrapolated Picard:
    cu     = 1.5*uh - 0.5*uh_old;
    cv     = 1.5*vh - 0.5*vh_old;
end



if (nonlinear_build_matrix==0)
      diag1 = [Omu;Omv];
      CD    = spdiags(diag1,0,Nu+Nv,Nu+Nv);
end


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

  
           
%% rhs 1


% convective terms, u-component
uIux       = Iu_ux*cu + yIu_ux; % convective velocity, u_bar
uAux       = Au_ux*uhn + yAu_ux;
yConv_ux   = Cux*( uIux.*uAux );

vIuy       = Iv_uy*cv + yIv_uy; % convective velocity, v_bar
uAuy       = Au_uy*uhn + yAu_uy;
yConv_uy   = Cuy*( vIuy.*uAuy );

% convective terms, v-component
uIvx       = Iu_vx*cu + yIu_vx; % convective velocity, u_bar  
vAvx       = Av_vx*vhn + yAv_vx;
yConv_vx   = Cvx*( uIvx.*vAvx );

vIvy       = Iv_vy*cv + yIv_vy; % convective velocity, v_bar
vAvy       = Av_vy*vhn + yAv_vy;   
yConv_vy   = Cvy*( vIvy.*vAvy );

% assemble rhs
f1(1:Nu)   = yConv_ux + yConv_uy + ...
             -Diffu*uhn - yDiffu - Fx + y_px + Gx*pn;

f1(Nu+1:Nu+Nv)  = yConv_vx + yConv_vy + ...
                -Diffv*vhn - yDiffv - Fy + y_py + Gy*pn;


C1         = Cux*spdiags(uIux,0,N1,N1);   
Conv_ux_11 = C1*Au_ux;

C1         = Cuy*spdiags(vIuy,0,N2,N2);   
Conv_uy_11 = C1*Au_uy;

C1         = Cvx*spdiags(uIvx,0,N3,N3);
Conv_vx_22 = C1*Av_vx;

C1         = Cvy*spdiags(vIvy,0,N4,N4);
Conv_vy_22 = C1*Av_vy;


% assemble matrix

CD_11      = - Diffu + Conv_ux_11 + Conv_uy_11;
CD_12      = spalloc(Nu,Nv,0);
CD_21      = spalloc(Nv,Nu,0);
CD_22      = - Diffv + Conv_vx_22 + Conv_vy_22; 

CD_1       = [CD_11 CD_12; CD_21 CD_22];

diag1      = [Omu;Omv];
CD         = spdiags(diag1,0,Nu+Nv,Nu+Nv);
CD         = CD + dt*a11*CD_1; 

      
% total matrix
Z         = [CD a11*dt*G; a11*dt*M Z2];
% also possible is:
% Z         = [CD a11*dt*G; a11*dt*M Z2];
% then we solve for kp instead of p; should not influence results for u


%% assemble total right-hand side
fM1       = M*V + yM;
f         = [f1;fM1];

% solve
q        = -Z\f;

ku1      = q(1:Nu);
kv1      = q(Nu+1:Nu+Nv);
kp       = q(Nu+Nv+1:Nu+Nv+Np);

% u1        = uhn + dt*a11*ku1;
% v1        = vhn + dt*a11*kv1;

uh        = uhn + dt*b1*ku1;
vh        = vhn + dt*b1*kv1;
p         = pn  + dt*b1*kp;
V         = [uh;vh];


if (BC_unsteady == 1)
% make V satisfy the incompressibility constraint at n+1 with a
% projection step
    t = tn + dt;
    boundary_conditions;
    interpolate_bc;
    operator_bc_divergence;
    t = tn;

    f       = M*V + yM;

    pressure_poisson;

    V       = V - Om_inv.*(G*dp);
    uh      = V(1:Nu);
    vh      = V(Nu+1:Nu+Nv);


end

if (EP==1)
    uh_old    = uhn;
    vh_old    = vhn;
end

% p_temp  = a_inv_11*(p - pn);
% p       = pn + b1*p_temp;

pressure_additional_solve;