function solver_unsteady_ROM_OpInf(snapshot_datas,dt_sample,t_sample,options)
%  Main solver file for unsteady calculations with reduced order model
%  based on operator inference

addpath('unsteady/ROM_bases_setup/');
addpath('opinf_stuff/');

%% load snapshot data for ROM basis construction
V_snapshots = [];

for i = 1:numel(snapshot_datas)
    snapshot_data = snapshot_datas(i);
    V_snapshots_ = load_snapshot_data("results/"+snapshot_data,dt_sample,t_sample);
    V_snapshots = [V_snapshots V_snapshots_];
end

%% construct ROM basis/bases
[basis,S] = Omega_POD(V_snapshots,options.rom.M,options);
options.rom.B = basis;

% semilogy(diag(S)/S(1,1))
% title("velocity singular values")

%% load snapshot data for operator inference
% reloading data could be avoided by storing a marker indicating which
% snapshot belongs to which trajectory -> information required because time
% difference quotients should only be taken within trajectories

A_dots = [];
As = [];
for i = 1:numel(snapshot_datas)
    snapshot_data = snapshot_datas(i);
    V_snapshots_ = load_snapshot_data("results/"+snapshot_data,dt_sample,t_sample);

    A_raw = options.rom.B'*(options.grid.Om.*V_snapshots_);
    [A_dot,A] = time_difference_quotient(A_raw, "forward euler",options); 
    A_dots = [A_dots A_dot];
    As = [As A];
end

options.rom.A = As;
options.rom.A_dot = A_dots;

%% construct ROM operators
% project velocity snapshots onto POD basis to get ROM coefficient
% snapshots for operator inference

% options = operator_rom(options);

[Diff_intrusive,Conv_intrusive] = rom_operator_wrapper(options,"intrusive");
[Diff_OpInf,Conv_OpInf] = rom_operator_wrapper(options,"OpInf");

Conv_intrusive_r = reduced_convection_operator(Conv_intrusive);

norm(Diff_intrusive - Diff_OpInf)
norm(reduced_convection_operator(Conv_intrusive)-reduced_convection_operator(-Conv_OpInf))

r = 6

%% initialize ROM solution

%% time stepping
