function solver_unsteady_ROM_OpInf(snapshot_datas,dt_sample,t_sample,options)
%  Main solver file for unsteady calculations with reduced order model
%  based on operator inference

addpath('unsteady/ROM_bases_setup/');
addpath('opinf_stuff/');

M = options.rom.M;
%% things that at best should be in parameter file

% Ms = [20:10:80];

Ms = M;
% Ms = 1:M;
% Ms = 10:10:M;

NMs = max(Ms);

% skip = 50;
skip = 1;


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

% [Diff_intrusive_,Conv_intrusive_] = rom_operator_wrapper(options,"intrusive");
% [Diff_intrusive_2,Conv_intrusive_2] = rom_operator_wrapper(options,"intrusive+");

[Diff_intrusive_,Conv_intrusive_] = rom_operator_wrapper(options,"intrusive+");

Conv_intrusive_tensor = reshape(Conv_intrusive_,M,M,M);
% Conv_intrusive_tensor2 = reshape(Conv_intrusive_2,M,M,M);

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
    % skip = 1;
    [A_dot,A] = time_difference_quotient(A_raw, "forward euler",options.time.dt,skip); 
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


% diff_errors = zeros(NMs,1);
% conv_errors = zeros(NMs,1);
% 
% diff_errors_ROM = zeros(NMs,1);
% conv_errors_ROM = zeros(NMs,1);
% 
% avg_rel_state_errors = zeros(NMs,1);
% avg_rel_state_errors_ROM = zeros(NMs,1);
% avg_rel_state_errors_intrusive = zeros(NMs,1);
% 
% three_term_errors = zeros(NMs,1);
% three_term_errors_ROM = zeros(NMs,1);
% three_term_errors_intrusive = zeros(NMs,1);
% 
% block_skewsymm_errors = zeros(NMs,1);
% block_skewsymm_errors_ROM = zeros(NMs,1);
% block_skewsymm_errors_intrusive = zeros(NMs,1);
% 
% rel_state_errors_intrusive2 = zeros(NMs,1);
% three_term_errors_intrusive2 = zeros(NMs,1);
% block_skewsymm_errors_intrusive2 = zeros(NMs,1);

rel_diff_errors = zeros(NMs,3,2); % [standard opinf, perm opinf, skew opinf] x [FOM projection data, closure-clean data]
rel_conv_errors = zeros(NMs,3,2);
rel_red_conv_errors = zeros(NMs,3,2);

avg_rel_state_error = zeros(NMs,3,2);

zero_perm_sum = zeros(NMs,3,2);
block_skewsymm = zeros(NMs,3,2);

avg_rel_state_errors_intrusive = zeros(NMs,1);
avg_rel_state_errors_best_approx = zeros(NMs,1);

zero_perm_sum_intrusive = zeros(NMs,1);
block_skewsymm_intrusive = zeros(NMs,1);

%% just plotting

% [~,S] = svd([As; vectorwise_kron(As)]);
%         figure(111)
%         s = diag(S);
%         plot(s/s(1),"kd-","DisplayName","FOM projection data")
%         hold on
%         sings_intrusive = s(1:100);
% 
%         sings_s = zeros(100,NMs);

%% construct ROM operators non-intrusively
a0 = A_raw(:,1);

    % rom_types = ["OpInf", "EC-OpInf Koike", "EC-OpInf skew"];
    rom_types = [ "EC-OpInf skew"];
    % labels = ["stan" "perm" "skew"];
    labels = ["skew"];
    n_rom_types = numel(rom_types);

    rel_diff_errors = zeros(NMs,n_rom_types,2); % [standard opinf, perm opinf, skew opinf] x [FOM projection data, closure-clean data]
    rel_conv_errors = zeros(NMs,n_rom_types,2);
    rel_red_conv_errors = zeros(NMs,n_rom_types,2);

    avg_rel_state_errors = zeros(NMs,n_rom_types,2);

    zero_perm_sum = zeros(NMs,n_rom_types,2);
    block_skewsymm = zeros(NMs,n_rom_types,2);


    zero_perm_sum_intrusive = zeros(NMs,1);

    block_skewsymm_intrusive = zeros(NMs,1);

for M_ = Ms

    options.rom.M = M_;
    options.rom.B = basis(:,1:M_);
    options.rom.A = As(1:M_,:);
    options.rom.A_dot = A_dots(1:M_,:);

    %% intrusive+ : intrusive with actually skew-symmetric convection operator
    % Diff_intrusive2 = Diff_intrusive_2(1:M_,1:M_);
    % Conv_intrusive2 = Conv_intrusive_tensor2(1:M_,1:M_,1:M_);
    % Conv_intrusive2 = Conv_intrusive2(:,:);
    %
    % A_intrusive2 = ROM_sim(Diff_intrusive2, -Conv_intrusive2,a0(1:M_),options.time.dt,size(A_raw,2));
    % rel_state_errors_intrusive2(M_) = relative_state_error(V_snapshots_,A_intrusive2,basis(:,1:M_));
    % three_term_constraint_ = three_term_prop_constraint(M_);
    % three_term_errors_intrusive2(M_) = norm(three_term_constraint_*reshape(Conv_intrusive2',M_^3,1));
    % block_skewsymm_constraint_ = block_skewsymm_constraint(M_);
    % block_skewsymm_errors_intrusive2(M_) = norm(block_skewsymm_constraint_*reshape(Conv_intrusive2',M_^3,1));

    %%

    Diff_intrusive = Diff_intrusive_(1:M_,1:M_);
    Conv_intrusive = Conv_intrusive_tensor(1:M_,1:M_,1:M_);
    Conv_intrusive = Conv_intrusive(:,:);

    % Conv_intrusive_r = reduced_convection_operator(Conv_intrusive);
    % Conv_OpInf_r = reduced_convection_operator(Conv_OpInf);

    % diff_error = norm(Diff_intrusive - Diff_OpInf)/norm(Diff_intrusive);
    % conv_error = norm(reduced_convection_operator(Conv_intrusive)-reduced_convection_operator(-Conv_OpInf))/norm(reduced_convection_operator(Conv_intrusive));
    % rel_diff_error = @(Diff_intrusive,Diff) norm(Diff_intrusive-Diff)/norm(Diff_intrusive);
    % rel_conv_error = @(Conv_intrusive,Conv) norm(Conv_intrusive-Conv)/norm(Conv_intrusive);
    % rel_red_conv_error = @(Conv_intrusive,Conv) norm(reduced_convection_operator(Conv_intrusive)-reduced_convection_operator(Conv))/norm(reduced_convection_operator(Conv_intrusive));
    rel_diff_error = @(Diff_intrusive,Diff) norm(Diff_intrusive-Diff,"fro")/norm(Diff_intrusive,"fro");
    rel_conv_error = @(Conv_intrusive,Conv) norm(Conv_intrusive-Conv,"fro")/norm(Conv_intrusive,"fro");
    rel_red_conv_error = @(Conv_intrusive,Conv) norm(reduced_convection_operator(Conv_intrusive-Conv),"fro")/norm(reduced_convection_operator(Conv_intrusive),"fro");

    % FOM projection data opinf
    % [Diff_stan_fp,Conv_stan_fp] = rom_operator_wrapper(options,"OpInf");
    % [Diff_perm_fp,Conv_perm_fp] = rom_operator_wrapper(options,"EC-OpInf Koike");
    % [Diff_skew_fp,Conv_skew_fp] = rom_operator_wrapper(options,"EC-OpInf skew");

    Diffs = zeros(M_,M_,n_rom_types,2);
    Convs = zeros(M_,M_^2,n_rom_types,2);

    for i = 1:n_rom_types
        [Diffs(:,:,i,1),Convs(:,:,i,1)] = rom_operator_wrapper(options,rom_types(i)); % 1 = FOM projection data
    end

    % diff_errors(M_) = diff_error;
    % conv_errors(M_) = conv_error;

    % A_intrusive = ROM_sim(Diff_intrusive, -Conv_intrusive,a0(1:M_),options.time.dt,size(A_raw,2));
    % avg_rel_state_errors_intrusive(M_) = relative_state_error(V_snapshots_,A_intrusive,basis(:,1:M_));
    avg_rel_state_errors_intrusive(M_) = average_relative_state_error(Diff_intrusive, -Conv_intrusive,a0s(1:M_,:),options.time.dt,size(A_raw,2),V_snapshots,basis(:,1:M_),skip);
    avg_rel_state_errors_best_approx(M_) = average_relative_state_error_best_approx(size(A_raw,2),V_snapshots,basis(:,1:M_),options);

    % A_opinf = ROM_sim(Diff_OpInf, Conv_OpInf,a0(1:M_),options.time.dt,size(A_raw,2));
    % avg_rel_state_errors(M_) = relative_state_error(V_snapshots_,A_opinf,basis(:,1:M_));
    % avg_rel_state_errors(M_) = average_relative_state_error(Diff_OpInf, Conv_OpInf,a0s(1:M_,:),options.time.dt,size(A_raw,2),V_snapshots,basis(:,1:M_));

    % avg_rel_state_errors(M_,1,1) = average_relative_state_error(Diff_stan_fp,Conv_stan_fp,a0s(1:M_,:),options.time.dt,size(A_raw,2),V_snapshots,basis(:,1:M_));
    % avg_rel_state_errors(M_,2,1) = average_relative_state_error(Diff_perm_fp,Conv_perm_fp,a0s(1:M_,:),options.time.dt,size(A_raw,2),V_snapshots,basis(:,1:M_));
    % avg_rel_state_errors(M_,3,1) = average_relative_state_error(Diff_skew_fp,Conv_skew_fp,a0s(1:M_,:),options.time.dt,size(A_raw,2),V_snapshots,basis(:,1:M_));

    %% ... and with closure-clean ROM simulation data
    [A_ROMs,A_dot_ROMs] = ROM_sims(Diff_intrusive, -Conv_intrusive,a0s(1:M_,:),options.time.dt,len_trajes,skip);

    options.rom.A = A_ROMs(1:M_,:);
    options.rom.A_dot = A_dot_ROMs(1:M_,:);

    %% just plotting

    % [~,S] = svd([A_ROMs; vectorwise_kron(A_ROMs)]);
    %     figure(111)
    %     s = diag(S);
    %     plot(s/s(1),"x-","DisplayName","closure-clean data r = "+num2str(M_))
    % 
    %     sings_s(1:(M_+M_^2),M_) = s(1:(M_+M_^2));
    %     legend("show")
    %     title("normalized singular value decay")
    % 
    %     ylabel("normalized magnitude")
    %     xlabel("singular value index")

    %%

    % [Diff_OpInf_ROM,Conv_OpInf_ROM] = rom_operator_wrapper(options,options.rom.rom_type);

    % [Diff_stan_cc,Conv_stan_cc] = rom_operator_wrapper(options,"OpInf");
    % [Diff_perm_cc,Conv_perm_cc] = rom_operator_wrapper(options,"EC-OpInf Koike");
    % [Diff_skew_cc,Conv_skew_cc] = rom_operator_wrapper(options,"EC-OpInf skew");

    for i = 1:n_rom_types
        [Diffs(:,:,i,2),Convs(:,:,i,2)] = rom_operator_wrapper(options,rom_types(i)); % 2 = closure-clean data
    end

    % Conv_OpInf_r_ROM = reduced_convection_operator(Conv_OpInf_ROM);


    % diff_error_ROM = norm(Diff_intrusive - Diff_OpInf_ROM)/norm(Diff_intrusive);
    % conv_error_ROM = norm(reduced_convection_operator(Conv_intrusive)-reduced_convection_operator(-Conv_OpInf_ROM))/norm(reduced_convection_operator(Conv_intrusive));
    %
    % diff_errors_ROM(M_) = diff_error_ROM;
    % conv_errors_ROM(M_) = conv_error_ROM;

    % A_opinf_ROM = ROM_sim(Diff_OpInf_ROM, Conv_OpInf_ROM,a0(1:M_),options.time.dt,size(A_raw,2));
    % avg_rel_state_errors_ROM(M_) = relative_state_error(V_snapshots_,A_opinf_ROM,basis(:,1:M_));
    % avg_rel_state_errors_ROM(M_) = average_relative_state_error(Diff_OpInf_ROM, Conv_OpInf_ROM,a0s(1:M_,:),options.time.dt,size(A_raw,2),V_snapshots,basis(:,1:M_));
    %%

    three_term_constraint_ = three_term_prop_constraint(M_);
    % % three_term_errors(M_) = norm(three_term_constraint_*reshape(Conv_OpInf',M_^3,1));
    % % three_term_errors_ROM(M_) = norm(three_term_constraint_*reshape(Conv_OpInf_ROM',M_^3,1));
    three_term_errors_intrusive(M_) = norm(three_term_constraint_*reshape(Conv_intrusive',M_^3,1));
    % 
    zero_perm_sum_intrusive(M_) = norm(three_term_constraint_*reshape(Conv_intrusive',M_^3,1));
    % 
    block_skewsymm_constraint_ = block_skewsymm_constraint(M_);
    % % block_skewsymm_errors(M_) = norm(block_skewsymm_constraint_*reshape(Conv_OpInf',M_^3,1));
    % % block_skewsymm_errors_ROM(M_) = norm(block_skewsymm_constraint_*reshape(Conv_OpInf_ROM',M_^3,1));
    % % block_skewsymm_errors_intrusive(M_) = norm(block_skewsymm_constraint_*reshape(Conv_intrusive',M_^3,1));
    % 
    block_skewsymm_intrusive(M_) = norm(block_skewsymm_constraint_*reshape(Conv_intrusive',M_^3,1));


    labels2 = ["FOM proj" "closure-clean"];
    for j = 1:2
        for i = 1:n_rom_types
            Diff = Diffs(:,:,i,j);
            % Diff = reshape(Diff,M_,M_);
            Conv = Convs(:,:,i,j);
            % Conv = reshape(Conv,M_,M_^2);

            rel_diff_errors(M_,i,j) = rel_diff_error(Diff_intrusive,Diff);
            rel_conv_errors(M_,i,j) = rel_conv_error(Conv_intrusive,-Conv);
            rel_red_conv_errors(M_,i,j) = rel_red_conv_error(Conv_intrusive,-Conv);

            avg_rel_state_errors(M_,i,j) = average_relative_state_error(Diff,Conv,a0s(1:M_,:),options.time.dt,size(A_raw,2),V_snapshots,basis(:,1:M_),skip);

            zero_perm_sum(M_,i,j) = norm(three_term_constraint_*reshape(Conv',M_^3,1));
            block_skewsymm(M_,i,j) = norm(block_skewsymm_constraint_*reshape(Conv',M_^3,1));

            labels_combined(i,j) = labels(i)+" "+ labels2(j);
            %%

        end
    end

end


%% plotting
markers = ["o" "+" "*"; "x", "s", "d"]

for i =1:n_rom_types

    % figure("rel operator errors FOM proj")
    figure(1)
    % semilogy(rel_diff_errors(:,i,1),'s-','DisplayName',"diffusion "+labels_combined(i,1))
    hold on
    semilogy(rel_conv_errors(:,i,1),'x-','DisplayName',"convection "+labels_combined(i,1))
    semilogy(rel_red_conv_errors(:,i,1),'d-','DisplayName',"reduced convection "+labels_combined(i,1))

    title("relative operator errors - FOM projection data")
    % legend("diffusion", "convection")
    legend("show")
    % title("original FOM data")

    % figure("rel operator errors closure-clean")
    figure(2)
    % semilogy(rel_diff_errors(:,i,2),'s-','DisplayName',"diffusion "+labels_combined(i,2))
    hold on
    semilogy(rel_conv_errors(:,i,2),'x-','DisplayName',"convection "+labels_combined(i,2))
    semilogy(rel_red_conv_errors(:,i,2),'+-','DisplayName',"reduced convection "+labels_combined(i,2))

    title("relative operator errors -closure-clean data")
    % legend("diffusion", "convection")
    legend("show")
    % title("closure-clean ROM data")

    for j =1:2

        figure(3)
        semilogy(avg_rel_state_errors(:,i,j),'d-','DisplayName',labels_combined(i,j))
        hold on
        % semilogy(avg_rel_state_errors_ROM,'x-')
        % semilogy(avg_rel_state_errors_intrusive,'o-')
        % semilogy(rel_state_errors_intrusive2,'<-')
        % legend("original FOM data","closure-clean data","intrusive")
        % legend("original FOM data","closure-clean data","intrusive","intrusive+")
        % title("relative state error (last trajectory only)")
        % title("average relative state error")
        legend("show")
        ylabel("average relative state error")
        xlabel("number of ROM modes")

        figure(4)
        % semilogy(three_term_errors(:,i,j),'d-','DisplayName',labels_combined(i,j))
        semilogy(zero_perm_sum(:,i,j),'d-','DisplayName',labels_combined(i,j))
        hold on
        % semilogy(three_term_errors_ROM,'x-')
        % semilogy(three_term_errors_intrusive,'o-')
        % semilogy(three_term_errors_intrusive2,'<-')
        % title("three term property error")
        % legend("original FOM data","closure-clean data","intrusive")
        % legend("original FOM data","closure-clean data","intrusive","intrusive+")
        title("zero permutation sum error")
        legend("show")

        figure(5)
        semilogy(block_skewsymm(:,i,j),'d-','DisplayName',labels_combined(i,j))
        hold on
        % semilogy(block_skewsymm_errors_ROM,'x-')
        % semilogy(block_skewsymm_errors_intrusive,'o-')
        % semilogy(block_skewsymm_errors_intrusive2,'<-')
        title("block skew-symmetry error")
        % legend("original FOM data","closure-clean data","intrusive")
        % legend("original FOM data","closure-clean data","intrusive","intrusive+")
        legend("show")

    end
end

energy_conservation_sims
eccomas_velocity_plot

r = 6

%% initialize ROM solution

%% time stepping
