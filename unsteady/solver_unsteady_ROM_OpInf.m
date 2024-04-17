function solver_unsteady_ROM_OpInf(snapshot_data,dt_sample,t_sample,options)
%  Main solver file for unsteady calculations with reduced order model
%  based on operator inference

%% load snapshot data
V_snapshots = load_snapshot_data(snapshot_data,dt_sample,t_sample);

%% construct ROM basis/bases
[basis,S] = Omega_POD(snapshot_matrix,options.rom.M,options);
options.rom.B = basis;

%% construct ROM operators

%% initialize ROM solution

%% time stepping
