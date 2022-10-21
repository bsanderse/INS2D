function V = get_velocity_from_streamfunction(psi,t,options)
% velocity values from streamfunction

Nx = options.grid.Nx;
Ny = options.grid.Ny;

% du/dy, like Su_uy
% note that we use hy here, and not gyd; this is to make sure that
% W' is a null vector of M
diag1  = 1./options.grid.hy;
W1D    = spdiags([-diag1 diag1],[-1 0],Ny,Ny);
W1D(1,end) = -W1D(1,1);
% extend to 2D
Wu_uy  = kron(W1D,speye(Nx));

% dv/dx, like Sv_vx
% note that we use hx here, and not gxd; this is to make sure that
% W' is a null vector of M
diag1  = 1./options.grid.hx;
W1D    = spdiags([-diag1 diag1],[-1 0],Nx,Nx);
W1D(1,end) = -W1D(1,1);
% extend to 2D
Wv_vx  = kron(speye(Ny),W1D);

% curl operator that acts on velocity fields
W = [-Wu_uy Wv_vx];
%         M_Null = W';
% gradient operator on psi is then given by W'
V = W' * psi;