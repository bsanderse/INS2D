function solver_unsteady_ROM_OpInf(snapshot_datas,dt_sample,t_sample,options)
%  Main solver file for unsteady calculations with reduced order model
%  based on operator inference

addpath('unsteady/ROM_bases_setup/');
addpath('opinf_stuff/');

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
% project velocity snapshots onto POD basis to get ROM coefficient
% snapshots for operator inference
if options.rom.rom_type == "OpInf" || options.rom.rom_type == "EC-OpInf"
    options.rom.A = options.rom.B'*(options.grid.Om.*V_snapshots);
end

options = operator_rom(options);

[Diff_intrusive,Conv_intrusive] = rom_operator_wrapper(options,"intrusive");
[Diff_OpInf,Conv_OpInf] = rom_operator_wrapper(options,"OpInf");

norm(Diff_intrusive - Diff_OpInf)
norm(reduced_convection_operator(Conv_intrusive)-Conv_OpInf)


%% initialize ROM solution

%% time stepping
