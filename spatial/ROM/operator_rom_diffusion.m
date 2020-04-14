function [yDiff,Diff] = operator_rom_diffusion(P,options)
% precompute convective operators
% projection with generic matrix P, size M x NV
% for momentum equation, P is B' or B'*Om_inv
% for Poisson equation, P is Bp'*M*Om_inv

NV = options.grid.Nu+options.grid.Nv;
Z  = zeros(NV,1);

B   = options.rom.B;
M   = options.rom.M;
Vbc = options.rom.Vbc;

% P can project to M or to Mp modes 
M1  = size(P,1);   

%% diffusion:
% the ROM discretization will read Diff*R + yDiff, where Diff is the reduced
% diffusion operator
% we get this by multiplying B'*Diffusion(B*R), where a call to Diffusion(V) returns
% D*V + yD
% this also works for the case that we have
% B'*Diffusion(B*R + Vbc) = 
% B'*(D*(B*R+Vbc) + yD) = (B'*D*B)*R + B'*D*Vbc + B'*yD =
% (B'*D*B)*R + B'*yDiff, where yDiff = B'*(D*Vbc + yD)
Diff = zeros(M1,M);

% first evaluate diffusion with zero velocity (or Vbc field) to only get boundary
% condition contribution:
% diffusion(Vbc) returns D*Vbc + yD
[d2u_bc,d2v_bc] = diffusion(Vbc,0,options,0);
yDiff           = P*[d2u_bc;d2v_bc];

% now evaluate each column of the reduced matrix as B'*D*B, where the
% boundary conditions are subtracted in order to get only the terms
% that involve the matrix
% we evaluate diffusion(B) = D*B + yD,
% and in order to get B'*D*B, we take
% diffusion(B) - diffusion(0)
[d2u_0,d2v_0] = diffusion(Z,0,options,0);
yDiff0        = P*[d2u_0;d2v_0];

for i=1:M
    [d2u,d2v] = diffusion(B(:,i),0,options,0);
    Diff(:,i) = P*[d2u;d2v] - yDiff0;
end