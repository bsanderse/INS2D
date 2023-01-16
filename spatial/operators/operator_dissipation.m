function options = operator_dissipation(options)

% boundary conditions
BC = options.BC;

Nx = options.grid.Nx;
Ny = options.grid.Ny;
Npx = options.grid.Npx;
Npy = options.grid.Npy;

% number of interior points and boundary points
Nu     = options.grid.Nu;
Nux_in = options.grid.Nux_in;
Nux_b  = options.grid.Nux_b;
Nux_t  = options.grid.Nux_t;
Nuy_in = options.grid.Nuy_in;
Nuy_b  = options.grid.Nuy_b;
Nuy_t  = options.grid.Nuy_t;
Nv     = options.grid.Nv;
Nvx_in = options.grid.Nvx_in;
Nvx_b  = options.grid.Nvx_b;
Nvx_t  = options.grid.Nvx_t;
Nvy_in = options.grid.Nvy_in;
Nvy_b  = options.grid.Nvy_b;
Nvy_t  = options.grid.Nvy_t;


Omu = options.grid.Omu;
Omv = options.grid.Omv;
hx  = options.grid.hx;
hy  = options.grid.hy;



%% take weighted average to get contribution to a certain velocity
% point
weight = 1/2;

% average dudx^2, effectively at u point
% multiply by finite volume size
diag1       = weight*ones(Nux_t,1);
A1D         = spdiags([diag1 diag1],[0 1],Nux_t-2,Nux_t-1);

% need to correct these dissipation terms because on "aligned" boundaries (d2udx2 and
% d2vdy2) we do not get a proper summation by parts identity
switch BC.u.left
    case 'dir'
        % "left boundary"
        A1D(1,1) = 1;
        %                 dudx2_u(1:Nux_in:end) = dudx2_u(1:Nux_in:Nu) +  ...
        %                     0.5*hy.*(uh(1:Nux_in:Nu).^2)/hx(1);
end
switch BC.u.right
    case 'dir'
        % "right boundary"
        A1D(end,end) = 1;
        %                 dudx2_u(Nux_in:Nux_in:end) = dudx2_u(Nux_in:Nux_in:Nu) + ...
        %                     0.5*hy.*(uh(Nux_in:Nux_in:Nu).^2)/hx(end);
end
Aux_u       = spdiags(Omu,0,Nu,Nu)*kron(speye(Nuy_in),A1D);
% dudx2_u     = Aux_u*(dudx.^2);

% average dudy^2, effectively at u point
% multiply by finite volume size
diag1       = weight*ones(Nuy_t,1);
A1D         = spdiags([diag1 diag1],[0 1],Nuy_t-2,Nuy_t-1);
Auy_u       = spdiags(Omu,0,Nu,Nu)*kron(A1D,speye(Nux_in));
% dudy2_u     = Auy_u*(dudy.^2);

% average dvdx^2, effectively at v point
% multiply by finite volume size
diag1       = weight*ones(Nvx_t,1);
A1D         = spdiags([diag1 diag1],[0 1],Nvx_t-2,Nvx_t-1);
Avx_v       = spdiags(Omv,0,Nv,Nv)*kron(speye(Nvy_in),A1D);
% dvdx2_v     = Avx_v*(dvdx.^2);

% average dvdy^2, effectively at v point
% multiply by finite volume size
diag1       = weight*ones(Nvy_t,1);
A1D         = spdiags([diag1 diag1],[0 1],Nvy_t-2,Nvy_t-1);

switch BC.v.low
    case 'dir'
        % "bottom boundary"
        A1D(1,1) = 1;
        %                 dvdy2_v(1:Nvx_in) = dvdy2_v(1:Nvx_in) + ...
        %                     0.5*hx.*(vh(1:Nvx_in).^2)/hy(1);
end
switch BC.v.up
    case 'dir'
        % "top boundary"
        A1D(end,end) = 1;
        %                 dvdy2_v(Nv-Nvx_in+1:end) = dvdy2_v(Nv-Nvx_in+1:end) + ...
        %                     0.5*hx.*(vh(Nv-Nvx_in+1:end).^2)/hy(end);
end
Avy_v       = spdiags(Omv,0,Nv,Nv)*kron(A1D,speye(Nvx_in));
% dvdy2_v     = Avy_v*(dvdy.^2);


%% then use definition of local kinetic energy to get Phi
% since the local KE  is defined on a pressure volume, as sum of
% its neighbouring velocity values (square), we get an additional
% interpolation step

% we need boundary contributions corresponding to ub^2, but we ignore these for now
% average from velocity point to pressure point
diag1     = weight*ones(Nux_t,1);
A1D       = spdiags([diag1 diag1],[0 1],Npx,Npx+1);
Au_k_BC  = ...
    BC_general(Npx+1,Nux_in,Npx+1-Nux_in,BC.u.left,BC.u.right,hx(1),hx(end));
% extend to 2D
Au_k      = kron(speye(Ny),A1D*Au_k_BC.B1D);

% we need boundary contributions for vb^2, ignore these for now
% average from velocity point to pressure point
diag1     = weight*ones(Nvy_t,1);
A1D       = spdiags([diag1 diag1],[0 1],Npy,Npy+1);
Av_k_BC  = ...
    BC_general(Npy+1,Nvy_in,Npy+1-Nvy_in,BC.v.low,BC.v.up,hy(1),hy(end));
% extend to 2D
Av_k      = kron(A1D*Av_k_BC.B1D,speye(Nx));


options.discretization.Aux_u = Aux_u;
options.discretization.Auy_u = Auy_u;
options.discretization.Avx_v = Avx_v;
options.discretization.Avy_v = Avy_v;
options.discretization.Au_k = Au_k;
options.discretization.Av_k = Av_k;

