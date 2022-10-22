function options = operator_temperature(options)
% construct convection and diffusion operators

% boundary conditions
BC = options.BC;

Npx = options.grid.Nx;
Npy = options.grid.Ny;

% number of interior points and boundary points
% Nu     = options.grid.Nu;
Nux_in = options.grid.Nux_in;
% Nux_b  = options.grid.Nux_b;
% Nux_t  = options.grid.Nux_t;
% Nuy_in = options.grid.Nuy_in;
% Nuy_b  = options.grid.Nuy_b;
% Nuy_t  = options.grid.Nuy_t;
% Nv     = options.grid.Nv;
% Nvx_in = options.grid.Nvx_in;
% Nvx_b  = options.grid.Nvx_b;
% Nvx_t  = options.grid.Nvx_t;
Nvy_in = options.grid.Nvy_in;
% Nvy_b  = options.grid.Nvy_b;
% Nvy_t  = options.grid.Nvy_t;
% 
hx  = options.grid.hx;
hy  = options.grid.hy;
% hxi = options.grid.hxi;
% hyi = options.grid.hyi;
% hxd = options.grid.hxd;
% hyd = options.grid.hyd;
% gxi = options.grid.gxi;
% gyi = options.grid.gyi;
gxd = options.grid.gxd;
gyd = options.grid.gyd;
% 
% Buvy = options.grid.Buvy;
% Bvux = options.grid.Bvux;
mat_hx = spdiags(hx,0,Npx,Npx);
mat_hy = spdiags(hy,0,Npy,Npy);


%%%%%%%%%%%%%%
% x-direction

%% differencing matrix
diag1  = ones(Npx,1);
C1D    = spdiags([-diag1 diag1],[0 1],Npx,Npx+1);
CTx    = kron(speye(Npy),C1D);

%% same but with face area
diag1    = ones(Npx,1); % Nkx = Npx
D1D      = spdiags([-diag1 diag1],[0 1],Npx,Npx+1);
% No BC
DTx      = kron(mat_hy,D1D);

%% interpolating u
I1D    = speye(Npx+1);

% boundary conditions
Iu_Tx_BC = BC_general(Npx+1,Nux_in,Npx+1-Nux_in,BC.u.left,BC.u.right,hx(1),hx(end));
% uLe_i  = interp1(y,uLe,yp);
% uRi_i  = interp1(y,uRi,yp);
% ybc    = kron(uLe_i,Iu_Tx_BC.ybc1) + kron(uRi_i,Iu_Tx_BC.ybc2);
% yIu_Tx = kron(mat_hy,I1D*Iu_Tx_BC.Btemp)*ybc;
Iu_Tx_BC.Bbc = kron(mat_hy,I1D*Iu_Tx_BC.Btemp);
Iu_Tx  = kron(mat_hy,I1D*Iu_Tx_BC.B1D);


%% averaging T
diag2  = 0.5*ones(Npx+1,1);
A1D    = spdiags([diag2 diag2],[0 1],Npx+1,Npx+2);

% boundary conditions
% [B1DT, BtempT, ybcl, ybcr] = 
AT_Tx_BC = BC_general_stag(Npx+2,Npx,2,BC.T.left,BC.T.right,hx(1),hx(end));
% ybcT    = kron(TLe,ybcl) + kron(TRi,ybcr);
% yAT_Tx  = kron(speye(Npy),A1D*AT_Tx_BC.BtempT)*ybcT;
AT_Tx   = kron(speye(Npy),A1D*AT_Tx_BC.B1D);
AT_Tx_BC.Bbc = kron(speye(Npy),A1D*AT_Tx_BC.Btemp);

%% differencing from centers to faces
diag3    = 1./gxd;
S1D      = spdiags([-diag3 diag3],[0 1],Npx+1,Npx+2);

% re-use BC generated for averaging T
DiffTx     = DTx * kron(speye(Npy),S1D*AT_Tx_BC.B1D);
DiffTx_BC  = DTx * kron(speye(Npy),S1D*AT_Tx_BC.Btemp);
% ySTx     = kron(speye(Npy),S1D*BtempT)*ybc;


%%%%%%%%%%%%%%
% y-direction

%% differencing matrix
diag1  = ones(Npy,1);
C1D    = spdiags([-diag1 diag1],[0 1],Npy,Npy+1);
CTy    = kron(C1D,speye(Npx));

%% same but with face area
diag1    = ones(Npy,1); % Nky = Npy
D1D      = spdiags([-diag1 diag1],[0 1],Npy,Npy+1);
% No BC
DTy      = kron(D1D,mat_hx);

%% interpolating v
I1D    = speye(Npy+1);

% boundary conditions
Iv_Ty_BC = BC_general(Npy+1,Nvy_in,Npy+1-Nvy_in, ...
                                     BC.v.low,BC.v.up,hy(1),hy(end));
% vLo_i  = interp1(x,vLo,xp);
% vUp_i  = interp1(x,vUp,xp);
% ybc    = kron(ybcl,vLo_i) + kron(ybcu,vUp_i);
% Iv_Ty.Bbc = kron(I1D*Iv_Ty_BC.Btemp,mat_hx)*ybc
Iv_Ty_BC.Bbc = kron(I1D*Iv_Ty_BC.Btemp,mat_hx);
Iv_Ty  = kron(I1D*Iv_Ty_BC.B1D,mat_hx);


%% averaging T
diag2  = 0.5*ones(Npy+1,1);
A1D    = spdiags([diag2 diag2],[0 1],Npy+1,Npy+2);

% boundary conditions
AT_Ty_BC = BC_general_stag(Npy+2,Npy,2, ...
                                     BC.T.low,BC.T.up,hy(1),hy(end));
% ybcT   = kron(ybcl,TLo) + kron(ybcr,TUp);
% yAT_Ty = kron(A1D*AT_Ty_BC.BtempT,speye(Npx))*ybcT;
AT_Ty  = kron(A1D*AT_Ty_BC.B1D,speye(Npx));
AT_Ty_BC.Bbc = kron(A1D*AT_Ty_BC.Btemp,speye(Npx));


%% differencing from centers to faces
diag3    = 1./gyd;
S1D      = spdiags([-diag3 diag3],[0 1],Npy+1,Npy+2);

% re-use BC generated for averaging t
DiffTy      = DTy*kron(S1D*AT_Ty_BC.B1D,speye(Npx));
DiffTy_BC   = DTy*kron(S1D*AT_Ty_BC.Btemp,speye(Npx));
% ySTy     = kron(S1D*BtempT,speye(Npx))*ybcT;


%% assemble total operator
DiffT    = DiffTx + DiffTy;
% yDiffT   = DTx*ySTx + DTy*ySTy;

%% store in options structure
options.discretization.CTx   = CTx;
options.discretization.CTy   = CTy;

options.discretization.Iu_Tx   = Iu_Tx;
options.discretization.Iv_Ty   = Iv_Ty;
options.discretization.AT_Tx   = AT_Tx;
options.discretization.AT_Ty   = AT_Ty;
options.discretization.Iu_Tx_BC   = Iu_Tx_BC;
options.discretization.Iv_Ty_BC   = Iv_Ty_BC;
options.discretization.AT_Tx_BC   = AT_Tx_BC;
options.discretization.AT_Ty_BC   = AT_Ty_BC;

% options.discretization.STx_BC   = STx_BC;
% options.discretization.STy_BC   = STy_BC;
options.discretization.DiffT    = DiffT;
options.discretization.DiffTx_BC = DiffTx_BC;
options.discretization.DiffTy_BC = DiffTy_BC;



