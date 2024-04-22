function solver_unsteady_ROM_OpInf(snapshot_datas,dt_sample,t_sample,options)
%  Main solver file for unsteady calculations with reduced order model
%  based on operator inference

%% load snapshot data
V_snapshots = [];
for i = 1:numel(snapshot_datas)
    snapshot_data = snapshot_datas(i);
    V_snapshots = [V_snapshots load_snapshot_data("results/"+snapshot_data,dt_sample,t_sample)];
end

%% construct ROM basis/bases
[basis,S] = Omega_POD(V_snapshots,options.rom.M,options);
options.rom.B = basis;

%% construct ROM operators

%% initialize ROM solution

%% time stepping
