function [u,v,p,options] = shear_layer_ROM_IC_opinf_basis(t,options)

%% note: this basis construction could be outsourced and only be done once
r = options.rom.M;

basis0 = eye(r);
basis = basis0;
for s = 1:r
    offset = zeros(r,1);
    offset(s) = 1; 
    basis = [basis, basis0 + offset]; % note: this basis includes redundant kronecker ambiguity pairs
end

%%
basis_j = basis(:,options.simulation_nr);

snapshot_data = options.rom.snapshot_data;
dt_sample = options.rom.dt_sample;
t_sample = options.rom.t_sample;

solver_unsteady_ROM_basis_construction;

V = options.rom.B*basis_j;

Nu = options.grid.Nu;
Npx = options.grid.Npx;
Npy = options.grid.Npy;

u   = V(1:Nu);
v   = V(Nu+1:end);
p   = zeros(Npx,Npy);

end