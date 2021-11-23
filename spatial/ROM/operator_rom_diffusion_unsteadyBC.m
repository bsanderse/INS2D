function [yDiff0,Diff,DiffBC] = operator_rom_diffusion_unsteadyBC(P,options)

B   = options.rom.B;
Vbc0 = options.rom.Vbc0;

Diffu  = options.discretization.Diffu;
Diffv  = options.discretization.Diffv;

Diff_   = P*blkdiag(Diffu, Diffv)*B;
DiffBC_ = P*blkdiag(Diffu, Diffv)*Vbc0;

yDiffu1 = options.bc_options1.discretization.yDiffu;
yDiffv1 = options.bc_options1.discretization.yDiffv;
yDiff0_1 = P*[yDiffu1; yDiffv1];

yDiffu2 = options.bc_options2.discretization.yDiffu;
yDiffv2 = options.bc_options2.discretization.yDiffv;
yDiff0_2 = P*[yDiffu2; yDiffv2];

Diff   = Diff_;
DiffBC = DiffBC_;
yDiff0 = [yDiff0_1, yDiff0_2];


