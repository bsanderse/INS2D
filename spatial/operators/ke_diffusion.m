% nu_T    = C_mu*(Akx*k)^2/(Aex*e)
% Dkx*(Anu*nu+C_mu*(Akx*k)^2/(Aex*e)).*Skx*k

%%%%%%%%%%%%%%
% x-direction

%% differencing from faces to centers
diag1    = ones(Npx,1); % Nkx = Npx
D1D      = spdiags([-diag1 diag1],[0 1],Npx,Npx+1);
% No BC
Dkx      = kron(mat_hy,D1D);


%% averaging from centers to faces
diag2    = 0.5*ones(Npx+1,1);
A1D      = spdiags([diag2 diag2],[0 1],Npx+1,Npx+2);

% BCs for k
% Ak_kx is already constructed in ke_convection
[B1Dk, Btempk, ybcl, ybcr] = BC_general_stag(Npx+2,Npx,2,BC.k.left,BC.k.right,hx(1),hx(end));
ybck     = kron(kLe,ybcl) + kron(kRi,ybcr);
yAk_kx   = kron(speye(Npy),A1D*Btempk)*ybck;
Ak_kx    = kron(speye(Npy),A1D*B1Dk);

% BCs for e
[B1De, Btempe, ybcl, ybcr] = BC_general_stag(Npx+2,Npx,2,BC.e.left,BC.e.right,hx(1),hx(end));
ybce     = kron(eLe,ybcl) + kron(eRi,ybcr);
yAe_ex   = kron(speye(Npy),A1D*Btempe)*ybce;
Ae_ex    = kron(speye(Npy),A1D*B1De);


%% differencing from centers to faces
diag3    = 1./gxd;
S1D      = spdiags([-diag3 diag3],[0 1],Npx+1,Npx+2);

% re-use BC generated for averaging k
Skx      = kron(speye(Npy),S1D*B1Dk);
ySkx     = kron(speye(Npy),S1D*Btempk)*ybck;

% re-use BC generated for averaging e
Sex      = kron(speye(Npy),S1D*B1De);
ySex     = kron(speye(Npy),S1D*Btempe)*ybce;


%%%%%%%%%%%%%%
% y-direction


%% differencing from faces to centers
diag1    = ones(Npy,1); % Nky = Npy
D1D      = spdiags([-diag1 diag1],[0 1],Npy,Npy+1);
% No BC
Dky      = kron(D1D,mat_hx);


%% averaging
diag2    = 0.5*ones(Npy+1,1);
A1D      = spdiags([diag2 diag2],[0 1],Npy+1,Npy+2);

% BCs for k:
% k is already constructed in ke_convection
[B1Dk, Btempk, ybcl, ybcu] = BC_general_stag(Npy+2,Npy,2,BC.k.low,BC.k.up,hy(1),hy(end));
ybck     = kron(ybcl,kLo) + kron(ybcu,kUp);
yAk_ky   = kron(A1D*Btempk,speye(Npx))*ybck;
Ak_ky    = kron(A1D*B1Dk,speye(Npx));

% BCs for e:
[B1De, Btempe, ybcl, ybcu] = BC_general_stag(Npy+2,Npy,2,BC.e.low,BC.e.up,hy(1),hy(end));
ybce     = kron(ybcl,eLo) + kron(ybcu,eUp);
yAe_ey   = kron(A1D*Btempe,speye(Npx))*ybce;
Ae_ey    = kron(A1D*B1De,speye(Npx));


%% differencing from centers to faces
diag3    = 1./gyd;
S1D      = spdiags([-diag3 diag3],[0 1],Npy+1,Npy+2);

% re-use BC generated for averaging k
Sky      = kron(S1D*B1Dk,speye(Npx));
ySky     = kron(S1D*Btempk,speye(Npx))*ybck;

% re-use BC generated for averaging e
Sey      = kron(S1D*B1De,speye(Npx));
ySey     = kron(S1D*Btempe,speye(Npx))*ybce;