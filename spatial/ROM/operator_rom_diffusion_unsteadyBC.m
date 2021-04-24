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
Diff = zeros(M1,M);
DiffBC = zeros(M1,Mbc);

[d2u_0,d2v_0] = mydiffusion(Z,0,options,0);
yDiff0        = P*[d2u_0;d2v_0];

for i=1:M
    [d2u,d2v] = mydiffusion(B(:,i),0,options,0);
    Diff(:,i) = P*[d2u;d2v] - yDiff0;
end

for i=1:Mbc
    [d2u,d2v] = mydiffusion(Vbc0(:,i),0,options,0);
    DiffBC(:,i) = P*[d2u;d2v] - yDiff0;
end

%% alternative

Diffu  = options.discretization.Diffu;
Diffv  = options.discretization.Diffv;
yDiffu = options.discretization.yDiffu;
yDiffv = options.discretization.yDiffv;

Diff_   = P*blkdiag(Diffu, Diffv)*B;
DiffBC_ = P*blkdiag(Diffu, Diffv)*Vbc0;
yDiff0_ = P*[yDiffu; yDiffv];

norm(Diff-Diff_)
norm(DiffBC-DiffBC_)
norm(yDiff0 - yDiff0_)

Diff   = Diff_;
DiffBC = DiffBC_;
yDiff0 = yDiffy0_;


