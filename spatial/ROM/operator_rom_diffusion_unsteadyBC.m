function [yDiff0,Diff,DiffBC] = operator_rom_diffusion_unsteadyBC(P,options)
% precompute convective operators
% projection with generic matrix P, size M x NV
% for momentum equation, P is B' or B'*Om_inv
% for Poisson equation, P is Bp'*M*Om_inv

NV = options.grid.Nu+options.grid.Nv;
Z  = zeros(NV,1);

B   = options.rom.B;
M   = options.rom.M;
Vbc0 = options.rom.Vbc0;
Mbc = size(Vbc0,2);

% P can project to M or to Mp modes
M1  = size(P,1);

% Phi^T(D_h(Phi a + Phi_bc a_bc) + y_D) = Diff a + DiffBC a_bc +
% y_Diff0

% Diff = zeros(M1,M);
% DiffBC = zeros(M1,Mbc);
% 
% [d2u_0,d2v_0] = mydiffusion(Z,0,options,0);
% yDiff0        = P*[d2u_0;d2v_0];
% 
% for i=1:M
%     [d2u,d2v] = mydiffusion(B(:,i),0,options,0);
%     Diff(:,i) = P*[d2u;d2v] - yDiff0;
% end
% 
% for i=1:Mbc
%     [d2u,d2v] = mydiffusion(Vbc0(:,i),0,options,0);
%     DiffBC(:,i) = P*[d2u;d2v] - yDiff0;
% end

%% alternative

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

yDiffu = options.discretization.yDiffu;
yDiffv = options.discretization.yDiffv;
yDiff0_ = P*[yDiffu; yDiffv];

norm(yDiff0_ - (yDiff0_1*options.rom.abc1(0)+yDiff0_2*options.rom.abc2(0)))

% norm(Diff-Diff_)
% norm(DiffBC-DiffBC_)
% norm(yDiff0 - (yDiff0_1*options.rom.abc1(0)+yDiff0_2*options.rom.abc2(0)))

Diff   = Diff_;
DiffBC = DiffBC_;
yDiff0 = [yDiff0_1, yDiff0_2];


