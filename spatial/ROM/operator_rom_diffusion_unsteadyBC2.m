function [yDiff_bc,Diff,Diff_inhom] = operator_rom_diffusion_unsteadyBC2(P,options)

B   = options.rom.B;

Diffu  = options.discretization.Diffu;
Diffv  = options.discretization.Diffv;

Diff_   = P*blkdiag(Diffu, Diffv)*B;

% if options.rom.rom_bc == 2 && options.rom.bc_recon == 3
if options.rom.bc_recon == 3 || options.rom.bc_recon == 5
    phi_bc = options.rom.phi_bc;
    M_bc = size(phi_bc,2);
    
    for i = 1:M_bc
        yBC = phi_bc(:,i);
        options = set_bc_vectors_from_yBC(options,yBC);
        
        yDiffu1 = options.discretization.yDiffu;
        yDiffv1 = options.discretization.yDiffv;
        yDiff_BC_(:,i) = P*[yDiffu1; yDiffv1];
    end
else
    yDiff_BC_ = -666;
end
   
if options.rom.bc_recon == 3
    phi_inhom = options.rom.phi_inhom;
    
    Diff_inhom_ = P*blkdiag(Diffu, Diffv)*phi_inhom;
else
    Diff_inhom_ = -666;
end

Diff   = Diff_;
Diff_inhom = Diff_inhom_;
yDiff_bc = yDiff_BC_;


