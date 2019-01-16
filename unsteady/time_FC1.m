% fully (energy) conservative time stepping
% solution of linear (!) saddlepoint system
% not reversible in time

CN     = 0.5;
f      = zeros(Nu+Nv+Np,1);
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
% also possible for the pressure?

% evaluate BC and force at intermediate time 
% should we also take the boundary conditions as extrapolations?
t = tn + 0.5*dt;
boundary_conditions;
interpolate_bc;
operator_bc_momentum;
force;

% evaluate BC for divergence at next time level
t = tn + dt;
boundary_conditions;
interpolate_bc;
operator_bc_divergence;

t = tn;

%% convective terms, u-component
% c^n * u^(n+1), c=u
uIux       = Iu_ux*cu+yIu_ux;                     % convective velocity, u_bar
uAux       = Au_ux*uhn+yAu_ux;
C1         = Cux*spdiags(uIux,0,N1,N1);   
Conv_ux_11 = C1*Au_ux;
yConv_ux   = Cux*( uIux.*uAux );

vIuy       = Iv_uy*cv+yIv_uy;                     % convective velocity, v_bar
uAuy       = Au_uy*uhn+yAu_uy;
C1         = Cuy*spdiags(vIuy,0,N2,N2);   
Conv_uy_11 = C1*Au_uy;
yConv_uy   = Cuy*( vIuy.*uAuy );


%% convective terms, v-component
uIvx       = Iu_vx*cu+yIu_vx;                 % convective velocity, u_bar  
vAvx       = Av_vx*vhn+yAv_vx;
C1         = Cvx*spdiags(uIvx,0,N3,N3);
Conv_vx_22 = C1*Av_vx;
yConv_vx   = Cvx*( uIvx.*vAvx );

vIvy       = Iv_vy*cv+yIv_vy;                 % convective velocity, v_bar
vAvy       = Av_vy*vhn+yAv_vy;
C1         = Cvy*spdiags(vIvy,0,N4,N4);
Conv_vy_22 = C1*Av_vy;
yConv_vy   = Cvy*( vIvy.*vAvy );


%% construct matrix (saddlepoint structure)
CD_11      = spdiags(Omu,0,Nu,Nu)/dt + ...
           CN*(- Diffu + Conv_ux_11 + Conv_uy_11);
CD_12      = spalloc(Nu,Nv,0);
CD_21      = spalloc(Nv,Nu,0);
CD_22      = spdiags(Omv,0,Nv,Nv)/dt + ...
           CN*(- Diffv + Conv_vx_22 + Conv_vy_22);    
CD         = [CD_11 CD_12; CD_21 CD_22];
Z          = [CD CN*G; M Z2]; 
% we can also multiply M by CN, then this should also be done with the residual M*V+yM


% right-hand side, contains residual; note there is no Omega*uh/dt involved
% because we solve for du
% are we missing boundary conditions for du here, which can be nonzero for
% unsteady bc?
f(1:Nu)          = yConv_ux + yConv_uy + ... 
                 -Diffu*uhn - yDiffu - Fx + Gx*p + y_px;
f(Nu+1:Nu+Nv)    = yConv_vx + yConv_vy  + ...
                 -Diffv*vhn - yDiffv - Fy + Gy*p + y_py;  
f(Nu+Nv+1:end)   = M*V + yM;

% solve
dq        = -Z\f;

dV        = dq(1:Nu+Nv);
dp        = dq(Nu+Nv+1:end);
du        = dV(1:Nu);
dv        = dV(Nu+1:end);

V         = V + dV;
uh        = uh + du;
vh        = vh + dv;
p         = p + dp;

if (EP==1)
    uh_old    = uhn;
    vh_old    = vhn;
end