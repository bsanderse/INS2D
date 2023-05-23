function yBC = get_bc_vector_yBC(t,options)
% construct boundary conditions

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
xp  = options.grid.xp;
yp  = options.grid.yp;


%% get BC values
uLo      = uBC(x,y(1),t,options);        % Nx
uUp      = uBC(x,y(end),t,options);      % Nx
% uLe      = uBC(x(1),y,t,options);
% uRi      = uBC(x(end),y,t,options);

uLo_i    = uBC(xin,y(1),t,options);      % Nxin
uUp_i    = uBC(xin,y(end),t,options);    % Nxin
uLe_i    = uBC(x(1),yp,t,options);       % Nyp
uRi_i    = uBC(x(end),yp,t,options);     % Nyp

% vLo      = vBC(x,y(1),t,options);
% vUp      = vBC(x,y(end),t,options);
vLe      = vBC(x(1),y,t,options);        % Ny
vRi      = vBC(x(end),y,t,options);      % Ny

vLo_i    = vBC(xp,y(1),t,options);       % Nxp
vUp_i    = vBC(xp,y(end),t,options);     % Nxp
vLe_i    = vBC(x(1),yin,t,options);      % Nyin
vRi_i    = vBC(x(end),yin,t,options);    % Nyin

yBC = [uLo(:); uUp(:); uLo_i(:); uUp_i(:); uLe_i(:); uRi_i(:);
       vLe(:); vRi(:); vLo_i(:); vUp_i(:); vLe_i(:); vRi_i(:)];


