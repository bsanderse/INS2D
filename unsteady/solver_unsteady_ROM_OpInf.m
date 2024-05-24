function solver_unsteady_ROM_OpInf(snapshot_datas,dt_sample,t_sample,options)
%  Main solver file for unsteady calculations with reduced order model
%  based on operator inference

addpath('unsteady/ROM_bases_setup/');
addpath('opinf_stuff/');

%% load snapshot data for ROM basis construction
V_snapshots = [];

for i = 1:numel(snapshot_datas)
% for i = 1:1
    snapshot_data = snapshot_datas(i);
    V_snapshots_ = load_snapshot_data("results/"+snapshot_data,dt_sample,t_sample);
    V_snapshots = [V_snapshots V_snapshots_];
end

%% construct ROM basis/bases
[basis,S] = Omega_POD(V_snapshots,options.rom.M,options);
options.rom.B = basis;

% semilogy(diag(S)/S(1,1))
% title("velocity singular values")

%% construct ROM operators intrusively

[Diff_intrusive_,Conv_intrusive_] = rom_operator_wrapper(options,"intrusive");

M = options.rom.M;
Conv_intrusive_tensor = reshape(Conv_intrusive_,M,M,M);

%% load snapshot data for operator inference
% project velocity snapshots onto POD basis to get ROM coefficient
% snapshots for operator inference

% reloading data could be avoided by storing a marker indicating which
% snapshot belongs to which trajectory -> information required because time
% difference quotients should only be taken within trajectories

A_dots = [];
As = [];

A_dots_ROM = [];
As_ROM = [];

n_trajes = numel(snapshot_datas); % number of trajectories
len_trajes = size(V_snapshots_,2); % length of trajectories
a0s = zeros(M,n_trajes);

for i = 1:n_trajes
    snapshot_data = snapshot_datas(i);
    V_snapshots_ = load_snapshot_data("results/"+snapshot_data,dt_sample,t_sample);

    A_raw = options.rom.B'*(options.grid.Om.*V_snapshots_);
    [A_dot,A] = time_difference_quotient(A_raw, "forward euler",options.time.dt); 
    A_dots = [A_dots A_dot];
    As = [As A];
    
    % % alternatively to projecting FOM data, perform ROM simulation to
    % % obtain ROM coefficient snapshot data without closure error
    % A_raw_ROM = ROM_sim(Diff_intrusive_, -Conv_intrusive_,A_raw(:,1),options.time.dt,size(A_raw,2));
    % % A_diff = vecnorm(A_raw-A_raw_ROM);
    % 
    % [A_dot_ROM,A_ROM] = time_difference_quotient(A_raw_ROM, "forward euler",options.time.dt); 
    % A_dots_ROM = [A_dots_ROM A_dot_ROM];
    % As_ROM = [As_ROM A_ROM];

    a0s(:,i) = A_raw(:,1);
end

% [A_ROMs,A_dot_ROMs] = ROM_sims(Diff_intrusive_, -Conv_intrusive_,a0s,options.time.dt,len_trajes);
% norm(A_ROMs-As_ROM) %unit test
% norm(A_dot_ROMs - A_dots_ROM)% unit test


Ms = [20:10:80];
% Ms = M;
% Ms = 1:M;
NMs = max(Ms);

diff_errors = zeros(NMs,1);
conv_errors = zeros(NMs,1);

diff_errors_ROM = zeros(NMs,1);
conv_errors_ROM = zeros(NMs,1);

rel_state_errors = zeros(NMs,1);
rel_state_errors_ROM = zeros(NMs,1);
rel_state_errors_intrusive = zeros(NMs,1);

%% construct ROM operators non-intrusively
a0 = A_raw(:,1);

for M_ = Ms

options.rom.M = M_;
options.rom.B = basis(:,1:M_);
options.rom.A = As(1:M_,:);
options.rom.A_dot = A_dots(1:M_,:);

[Diff_OpInf,Conv_OpInf] = rom_operator_wrapper(options,"OpInf");

Diff_intrusive = Diff_intrusive_(1:M_,1:M_);
Conv_intrusive = Conv_intrusive_tensor(1:M_,1:M_,1:M_);
Conv_intrusive = Conv_intrusive(:,:);

Conv_intrusive_r = reduced_convection_operator(Conv_intrusive);
Conv_OpInf_r = reduced_convection_operator(Conv_OpInf);

diff_error = norm(Diff_intrusive - Diff_OpInf)/norm(Diff_intrusive);
conv_error = norm(reduced_convection_operator(Conv_intrusive)-reduced_convection_operator(-Conv_OpInf))/norm(reduced_convection_operator(Conv_intrusive));

diff_errors(M_) = diff_error;
conv_errors(M_) = conv_error;

A_intrusive = ROM_sim(Diff_intrusive, -Conv_intrusive,a0(1:M_),options.time.dt,size(A_raw,2));
rel_state_errors_intrusive(M_) = relative_state_error(V_snapshots_,A_intrusive,basis(:,1:M_));

A_opinf = ROM_sim(Diff_OpInf, Conv_OpInf,a0(1:M_),options.time.dt,size(A_raw,2));
rel_state_errors(M_) = relative_state_error(V_snapshots_,A_opinf,basis(:,1:M_));

%% ... and with closure-clean ROM simulation data
[A_ROMs,A_dot_ROMs] = ROM_sims(Diff_intrusive, -Conv_intrusive,a0s(1:M_,:),options.time.dt,len_trajes);

options.rom.A = A_ROMs(1:M_,:);
options.rom.A_dot = A_dot_ROMs(1:M_,:);

[Diff_OpInf_ROM,Conv_OpInf_ROM] = rom_operator_wrapper(options,"OpInf");

Conv_OpInf_r_ROM = reduced_convection_operator(Conv_OpInf_ROM);


diff_error_ROM = norm(Diff_intrusive - Diff_OpInf_ROM)/norm(Diff_intrusive);
conv_error_ROM = norm(reduced_convection_operator(Conv_intrusive)-reduced_convection_operator(-Conv_OpInf_ROM))/norm(reduced_convection_operator(Conv_intrusive));

diff_errors_ROM(M_) = diff_error_ROM;
conv_errors_ROM(M_) = conv_error_ROM;

A_opinf_ROM = ROM_sim(Diff_OpInf_ROM, Conv_OpInf_ROM,a0(1:M_),options.time.dt,size(A_raw,2));
rel_state_errors_ROM(M_) = relative_state_error(V_snapshots_,A_opinf_ROM,basis(:,1:M_));
%%


end

figure
semilogy(diff_errors,'d-')
hold on
semilogy(conv_errors,'x-')

title("relative operator errors")
legend("diffusion", "convection")
title("original FOM data")

figure
semilogy(diff_errors_ROM)
hold on
semilogy(conv_errors_ROM,'x-')

title("relative operator errors")
legend("diffusion", "convection")
title("closure-clean ROM data")

figure
semilogy(rel_state_errors,'d-')
hold on
semilogy(rel_state_errors_ROM,'x-')
semilogy(rel_state_errors_intrusive,'o-')
legend("original FOM data","closure-clean data","intrusive")
title("relative state error (last trajectory only)")

r = 6

%% initialize ROM solution

%% time stepping
