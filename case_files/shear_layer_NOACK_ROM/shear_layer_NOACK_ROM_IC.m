function [u,v,p,options] = shear_layer_NOACK_ROM_IC(t,options)
% initial velocity field Taylor-Green

Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;
Npx = options.grid.Npx;
Npy = options.grid.Npy;


xu = options.grid.xu;
xv = options.grid.xv;
yu = options.grid.yu;
yv = options.grid.yv;


delta = options.fluid.d_layer;

U1 = options.fluid.U1;
U2 = options.fluid.U2;
Uavg = (U1+U2)/2;

pert = 0.01*Uavg;

% base flow:
u = 0.5*(U1 + U2) + 0.5*(U1-U2)*tanh(yu./delta);
v = zeros(Nvx_in,Nvy_in);
p = zeros(Npx,Npy);

% perturbations from Orr-Sommerfeld analysis
I      = sqrt(-1);
alpha  = 0.4070-0.0845*I; % see paper Noack, JFM 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get Orr-Sommerfeld eigenmodes as chebfuns
disp('determine initial condition from Orr-Sommerfeld analysis...')
[ mode_u, mode_v, omega ] = shear_layer_OrrSommerfeld( alpha, t ,options );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% the chebfuns can be evaluated at y locations
% U_perturb2D = kron(mode_u(yp),ones(Nux_in,1));
% V_perturb2D = kron(mode_v(yin),ones(Nvx_in,1));

u_pert = real(mode_u(yu).*exp(I*(alpha*xu - omega*t)));
v_pert = real(mode_v(yv).*exp(I*(alpha*xv - omega*t)));

% [u_pert, v_pert, omega] = shear_layer_OrrSommerfeld( alpha, t ,options );

pert = pert/max(real(mode_v));

u = u + pert*u_pert;
v = v + pert*v_pert;

% store in options structure, can be used in postprocessing or in BC
options.fluid.mode_u = mode_u;
options.fluid.mode_v = mode_v;
options.fluid.omega  = omega;
options.fluid.alpha  = alpha;
options.fluid.pert   = pert;

end