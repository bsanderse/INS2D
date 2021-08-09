function yBC = get_bc_vector_yBC(t,options)
% construct boundary conditions


%% get settings from options structure

% steady
steady = options.case.steady;
% 4th order
order4 = options.discretization.order4;
% Reynolds number
Re = options.fluid.Re;
% boundary conditions
BC = options.BC;

global uBC vBC dudtBC dvdtBC;

% type of stress tensor
visc = options.case.visc;


if (order4 == 1)
    alfa   = options.discretization.alfa;
end

%% grid settings

% number of interior points and boundary points
Nux_in = options.grid.Nux_in;
Nvy_in = options.grid.Nvy_in;
Np     = options.grid.Np;
Npx    = options.grid.Npx;
Npy    = options.grid.Npy;


xin = options.grid.xin;
yin = options.grid.yin;
x   = options.grid.x;
y   = options.grid.y;
hx  = options.grid.hx;
hy  = options.grid.hy;
xp  = options.grid.xp;
yp  = options.grid.yp;


%% get BC values
uLo      = uBC(x,y(1),t,options);
uUp      = uBC(x,y(end),t,options);
% uLe      = uBC(x(1),y,t,options);
% uRi      = uBC(x(end),y,t,options);

uLo_i    = uBC(xin,y(1),t,options);
uUp_i    = uBC(xin,y(end),t,options);
uLe_i    = uBC(x(1),yp,t,options);
uRi_i    = uBC(x(end),yp,t,options);

% vLo      = vBC(x,y(1),t,options);
% vUp      = vBC(x,y(end),t,options);
vLe      = vBC(x(1),y,t,options);
vRi      = vBC(x(end),y,t,options);

vLo_i    = vBC(xp,y(1),t,options);
vUp_i    = vBC(xp,y(end),t,options);
vLe_i    = vBC(x(1),yin,t,options);
vRi_i    = vBC(x(end),yin,t,options);

yBC = [uLo(:); uUp(:); uLo_i(:); uUp_i(:); uLe_i(:); uRi_i(:);
       vLe(:); vRi(:); vLo_i(:); vUp_i(:); vLe_i(:); vRi_i(:)];


