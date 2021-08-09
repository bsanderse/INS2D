function [yDiff_bc,Diff,Diff_inhom] = operator_rom_diffusion_unsteadyBC2(P,options)
% precompute convective operators
% projection with generic matrix P, size M x NV
% for momentum equation, P is B' or B'*Om_inv
% for Poisson equation, P is Bp'*M*Om_inv

B   = options.rom.B;
phi_inhom = options.rom.phi_inhom;

Diffu  = options.discretization.Diffu;
Diffv  = options.discretization.Diffv;

Diff_   = P*blkdiag(Diffu, Diffv)*B;
Diff_inhom_ = P*blkdiag(Diffu, Diffv)*phi_inhom;

phi_bc = options.rom.phi_bc;
M_bc = size(phi_bc,2);

for i = 1:M_bc
    yBC = phi_bc(:,i);
    options = set_bc_vectors_from_yBC(options,yBC);
    
    yDiffu1 = options.discretization.yDiffu;
    yDiffv1 = options.discretization.yDiffv;
    yDiff_BC_(:,i) = P*[yDiffu1; yDiffv1];
end 

Diff   = Diff_;
Diff_inhom = Diff_inhom_;
yDiff_bc = yDiff_BC_;


