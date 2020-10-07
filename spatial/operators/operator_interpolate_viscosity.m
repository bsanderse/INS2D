function options = operator_interpolate_viscosity(options)

% average (turbulent) viscosity to cell faces: from nu at xp,yp to nu at
% ux, uy, vx, vy locations

% see also ke_viscosity.m

% averaging weight:
weight = 1/2;

BC = options.BC;

Npx = options.grid.Nx;
Npy = options.grid.Ny;

% number of interior points and boundary points
Nux_in = options.grid.Nux_in;
Nuy_in = options.grid.Nuy_in;
Nvx_in = options.grid.Nvx_in;
Nvy_in = options.grid.Nvy_in;


hx  = options.grid.hx;
hy  = options.grid.hy;

Buvy = options.grid.Buvy;
Bvux = options.grid.Bvux;
Bkux = options.grid.Bkux;
Bkvy = options.grid.Bkvy;


% set BC for nu
% in the periodic case, the value of nu is not needed
% in all other cases, homogeneous (zero) Neumann conditions are used

if (strcmp(BC.u.left,'per') && strcmp(BC.u.right,'per'))
    BC.nu.left  = 'per';
    BC.nu.right = 'per';
else
    BC.nu.left  = 'sym';
    BC.nu.right = 'sym';
end

if (strcmp(BC.v.low,'per') && strcmp(BC.v.up,'per'))
    BC.nu.low = 'per';
    BC.nu.up  = 'per';
else
    BC.nu.low = 'sym';
    BC.nu.up  = 'sym';
end



%% nu to ux positions

A1D = speye(Npx+2);
A1D = Bkux*A1D;

% boundary conditions for nu; mapping from Npx (k) points to Npx+2 points
Anu_ux_BC = BC_general_stag(Npx+2,Npx,2, ...
                                      BC.nu.left,BC.nu.right,hx(1),hx(end));
% then map back to Nux_in+1 points (ux-faces) with Bkux
                                  
% extend to 2D
Anu_ux  = kron(speye(Nuy_in),A1D*Anu_ux_BC.B1D);
% ybc     = kron(nuLe,ybcl) + kron(nuRi,ybcr);
% yAnu_ux = kron(speye(Nuy_in),Bkux*A1D*Anu_ux_BC.Btemp)*ybc;
Anu_ux_BC.Bbc = kron(speye(Nuy_in),A1D*Anu_ux_BC.Btemp);

    
% so nu at ux is given by:
% (Anu_ux * nu + yAnu_ux)


%% nu to uy positions

% average in x-direction
diag1  = weight*ones(Npx+1,1);
A1D    = spdiags([diag1 diag1],[0 1],Npx+1,Npx+2);
% then map to Nux_in points (like Iv_uy) with Bvux
A1D    = Bvux*A1D;

% calculate average in y-direction; no boundary conditions
diag1  = weight*ones(Npy+1,1);
A1Dy   = spdiags([diag1 diag1],[0 1],Npy+1,Npy+2);
A2Dy   = kron(A1Dy,speye(Nux_in));

% boundary conditions for nu in x-direction; 
% mapping from Npx (nu) points to Npx+2 points
Anu_uy_BC_lr = BC_general_stag(Npx+2,Npx,2, ...
                                      BC.nu.left,BC.nu.right,hx(1),hx(end));
% extend BC to 2D
A2D    = kron(speye(Npy+2),A1D*Anu_uy_BC_lr.B1D);
% nuLe_i = [nuLe(1);nuLe;nuLe(end)];
% nuRi_i = [nuRi(1);nuRi;nuRi(end)];
% ybc    = kron(nuLe_i,Anu_uy_BC_lr.ybc1)+ kron(nuRi_i,Anu_uy_BC_lr.ybc2);     
% yAnu_uy_lr = kron(speye(Npy+2),A1D*Anu_uy_BC_lr.Btemp)*ybc;            

% apply bc in y-direction
Anu_uy_BC_lu = BC_general_stag(Npy+2,Npy,2,...
                                           BC.nu.low,BC.nu.up,hy(1),hy(end));
                                       
% ybc    = kron(Anu_uy_BC_lu.Btemp*Anu_uy_BC_lu.ybc1,nuLo) + kron(Anu_uy_BC_lu.Btemp*Anu_uy_BC_lu.ybc2,nuUp);
% yAnu_uy_lu = A2D*ybc;

A2Dx   = A2D*kron(Anu_uy_BC_lu.B1D,speye(Npx));            
             
Anu_uy  = A2Dy*A2Dx;
% yAnu_uy = A2Dy*(yAnu_uy_lu + yAnu_uy_lr);

% NEW:
Anu_uy_BC_lr.B2D = A2Dy*kron(speye(Npy+2),A1D*Anu_uy_BC_lr.Btemp);
% ybc        = kron(nuLe_i,Anu_uy_BC_lr.ybc1)+ kron(nuRi_i,Anu_uy_BC_lr.ybc2);     
% yAnu_uy_lr = Anu_uy_BC_lr.B2D*ybc;

Anu_uy_BC_lu.B2D = A2Dy*A2D*kron(Anu_uy_BC_lu.Btemp,speye(Npx));
% ybc         = kron(Anu_uy_BC_lu.ybc1,nuLo) + kron(Anu_uy_BC_lu.ybc2,nuUp);
% yAnu_uy_lu  = Anu_uy_BC_lu.B2D*ybc;

% so nu at uy is given by:
% (Anu_uy * nu + yAnu_uy)


%% nu to vx positions
diag1     = weight*ones(Npy+1,1);
A1D       = spdiags([diag1 diag1],[0 1],Npy+1,Npy+2);
% map to Nvy_in points (like Iu_vx) with Buvy
A1D       = Buvy*A1D;

% calculate average in x-direction; no boundary conditions
diag1     = weight*ones(Npx+1,1);
A1Dx      = spdiags([diag1 diag1],[0 1],Npx+1,Npx+2);
A2Dx      = kron(speye(Nvy_in),A1Dx);


% boundary conditions for nu in y-direction; 
% mapping from Npy (nu) points to Npy+2 points
Anu_vx_BC_lu = BC_general_stag(Npy+2,Npy,2,...
                                           BC.nu.low,BC.nu.up,hy(1),hy(end));
% extend BC to 2D
A2D       = kron(A1D*Anu_vx_BC_lu.B1D,speye(Npx+2));


% apply boundary conditions also in x-direction:
Anu_vx_BC_lr = BC_general_stag(Npx+2,Npx,2,...
                                           BC.nu.left,BC.nu.right,hx(1),hx(end));
                                     
A2Dy      = A2D*kron(speye(Npy),Anu_vx_BC_lr.B1D);

Anu_vx     = A2Dx*A2Dy;

% OLD:
% nuLo_i     = [nuLo(1);nuLo;nuLo(end)];
% nuUp_i     = [nuUp(1);nuUp;nuUp(end)];
% ybc        = kron(Anu_vx_BC_lu.ybc1,nuLo_i) + kron(Anu_vx_BC_lu.ybc2,nuUp_i);
% yAnu_vx_lu = kron(A1D*Anu_vx_BC_lu.Btemp,speye(Npx+2))*ybc;
% yAnu_vx_lu1 = A2Dx*yAnu_vx_lu;
% 
% ybc       = kron(nuLe,Anu_vx_BC_lr.Btemp*Anu_vx_BC_lr.ybc1) + kron(nuRi,Anu_vx_BC_lr.Btemp*Anu_vx_BC_lr.ybc2);
% yAnu_vx_lr1 = A2Dx*A2D*ybc;


% NEW:
Anu_vx_BC_lu.B2D = A2Dx*kron(A1D*Anu_vx_BC_lu.Btemp,speye(Npx+2));
Anu_vx_BC_lr.B2D = A2Dx*A2D*kron(speye(Npy),Anu_vx_BC_lr.Btemp);

% % in y-direction
% ybc        = kron(Anu_vx_BC_lu.ybc1,nuLo_i) + kron(Anu_vx_BC_lu.ybc2,nuUp_i);
% yAnu_vx_lu2 = Anu_vx_lu.B2D*ybc;
% % in x-direction
% ybc        = kron(nuLe,Anu_vx_BC_lr.ybc1) + kron(nuRi,Anu_vx_BC_lr.ybc2);
% yAnu_vx_lr2 = Anu_vx_lr.B2D*ybc;

% so nu at uy is given by:
% (Anu_vx * nu + yAnu_vx)


%% nu to vy positions
A1D       = speye(Npy+2);
% then map back to Nvy_in+1 points (vy-faces) with Bkvy                                 
A1D       = Bkvy*A1D;

% boundary conditions for nu; mapping from Npy (nu) points to Npy+2 (vy faces) points
Anu_vy_BC = BC_general_stag(Npy+2,Npy,2, ...
                                      BC.nu.low,BC.nu.up,hy(1),hy(end));
                                
% extend to 2D
Anu_vy  = kron(A1D*Anu_vy_BC.B1D,speye(Nvx_in));
Anu_vy_BC.Bbc = kron(A1D*Anu_vy_BC.Btemp,speye(Nvx_in));

     
% so nu at vy is given by:
% (Anu_vy * k + yAnu_vy)

%%
options.discretization.Anu_ux    = Anu_ux;
options.discretization.Anu_ux_BC = Anu_ux_BC;

options.discretization.Anu_uy    = Anu_uy;
options.discretization.Anu_uy_BC_lr = Anu_uy_BC_lr;
options.discretization.Anu_uy_BC_lu = Anu_uy_BC_lu; 

options.discretization.Anu_vx    = Anu_vx;
options.discretization.Anu_vx_BC_lr = Anu_vx_BC_lr;
options.discretization.Anu_vx_BC_lu = Anu_vx_BC_lu;

options.discretization.Anu_vy    = Anu_vy;
options.discretization.Anu_vy_BC = Anu_vy_BC;
