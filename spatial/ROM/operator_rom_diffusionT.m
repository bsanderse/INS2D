function [yDiffT,DiffT] = operator_rom_diffusionT(PT,options)
% function [DiffT] = operator_rom_diffusionT(PT,options)
% precompute convective operators
% projection with generic matrix P, size M x NV
% for momentum equation, P is B' or B'*Om_inv
% for Poisson equation, P is Bp'*M*Om_inv

NT = options.grid.NT;
ZT  = zeros(NT,1);

BT   = options.rom.BT;
MT   = options.rom.MT;
% Tbc = options.rom.Tbc;

% P can project to M or to Mp modes 
MT1  = size(PT,1);   

%% diffusion:
% the ROM discretization will read Diff*R + yDiff, where Diff is the reduced
% diffusion operator
% we get this by multiplying B'*Diffusion(B*R), where a call to Diffusion(V) returns
% D*V + yD
% this also works for the case that we have
% B'*Diffusion(B*R + Vbc) = 
% B'*(D*(B*R+Vbc) + yD) = (B'*D*B)*R + B'*D*Vbc + B'*yD =
% (B'*D*B)*R + B'*yDiff, where yDiff = B'*(D*Vbc + yD)
DiffT = zeros(MT1,MT);

% first evaluate diffusion with zero velocity (velocity at t=0 or Vbc field) to only get boundary
% condition contribution:
% diffusion(Vbc) returns D*Vbc + yD

%% Newly added for temperature for non-homogenous boundary conditions
%It should be in this way
% [d2T_bc] = diffusion_temperature(Tbc,0,options,0);
% yDiffT           = PT*[d2T_bc];
% But this is simply manupulated right now as
NT=options.grid.NT;
yDiffT = zeros(NT,1);

%% now evaluate each column of the reduced matrix as B'*D*B, where the
% boundary conditions are subtracted in order to get only the terms
% that involve the matrix
% we evaluate diffusion(B) = D*B + yD,
% and in order to get B'*D*B, we take
% diffusion(B) - diffusion(0)
[d2T_0] = diffusion_temperature(ZT,0,options,0);
yDiffT0        = PT*[d2T_0];

%%For testing only 
% diffphi = zeros(NT,NT);
% %testing end

for i=1:MT
    [d2T] = diffusion_temperature(BT(:,i),0,options,0);
%     %testing only
%     diffphi(:,i)=d2T-d2T_0;
% %     testing end;  
    DiffT(:,i) = PT*[d2T] - yDiffT0;
end
% % Testinng again  remove this after testing
% DiffT=diffphi;
% diffphi=reshape(DiffT,NT*NT,1);
% save("ROM_diffphi.dat", "diffphi", "-ascii");
% disp("Hi")