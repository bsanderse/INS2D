function [B,B_inv,M] = getVelocityBasisFDG(options)
% Finite Difference Galerkin = reformulation in streamfunction form

%         % generate null vectors
%         Null_u = zeros(options.grid.Nu,1);
%         Null_u(2)    = 1;
%         Null_u(2+Nx) = -1;
%         Null_v = zeros(options.grid.Nv,1);
%         Null_v(1+Nx) = -1;
%         Null_v(2+Nx) = 1;

% use vorticity operator to construct null space
% note: need to get rid of extra periodic entries

%% for periodic BC, not covering entire mesh

Nu  = options.grid.Nu;
Nv  = options.grid.Nv;

Nx  = options.grid.Nx;
Ny  = options.grid.Ny;

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
% null space is then given by W'
B = W';

% is the divergence of B indeed zero?
div_basis = max(abs(options.discretization.M*B),[],1); %
% max over all columns:
maxdiv_basis = max(div_basis);
if (maxdiv_basis > 1e-14)
    error(['FDG basis not divergence free: ' num2str(maxdiv_basis) '\n']);
%     warning('Adding basis for pressure\n');
%     div_free = 0;
% else
%     div_free = 1;
end


% we can thus expand V, such that M*V=M*B*psi=0 for any psi:
% V = B*psi
% and
% psi = (B'*Om*B)^{-1} * B^T * Om * V

Om     = options.grid.Om;
Om_mat = spdiags(Om,0,Nu+Nv,Nu+Nv);

% ROM dimension (=NV - Np)
% there is no real dimension reduction as in the case of
% truncation, it's just going from velocity to streamfunction space
M = size(B,2);


% get the oblique projection (since B is not
% orthonormal)
B_inv = decomposition(B.'*Om_mat*B);


% options.rom.B = B;
% options.rom.M = M;

% Vbc = zeros(Nu+Nv,1);
% options.rom.Vbc = Vbc;