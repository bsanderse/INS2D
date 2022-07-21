%  Main solver file for unsteady calculations with reduced order model

addpath('unsteady/ROM_bases_setup/');

%% load snapshot data
if (j==1) || changing_snapshotdata
    disp('loading data...');
    snapshots = load(snapshot_data,'uh_total','vh_total','p_total','dt','t_end','Re','k','umom','vmom','maxdiv','Vbc','options', ...
        'k_diff', 'k_conv', 'k_pres', 'k_force', 'k_diffBC', 'k_presBC','k_obc');
    % snapshots.U = [snapshots.uh_total; snapshots.vh_total];
    
    % dt that was used for creating the snapshot matrix:
    dt_snapshots = snapshots.dt;
    % velocity field as snapshot matrix:
    V_total_snapshots = [snapshots.uh_total';snapshots.vh_total'];
    
    pro_X_h = V_total_snapshots;
    
    % check input dimensions
    Nspace  = size(V_total_snapshots,1); % total number of unknowns (Nu+Nv) of the original model
    Nu      = options.grid.Nu;
    Nv      = options.grid.Nv;
    
    some_snapshot_checks;
    
    FOM_data_consistency_check(options,snapshots.options)
    
    % subtract non-homogeneous BC contribution
    inhom_substraction;
    
    if options.rom.rom_bc == 2
        X_inhom = Vbc;
    end
     
    snapshot_sample = get_snapshot_sample(dt_sample,dt_snapshots,t_sample);
    
    % select snapshots
    V_svd = V_total_snapshots(:,snapshot_sample);
    
    X_h = pro_X_h(:,snapshot_sample);
    
    clear V_total_snapshots;  
    
    X_bc = get_X_bc(t_start,t_end,dt,snapshot_sample,options);
end

options.rom.bases_construction = "mthesis";
% options.rom.bases_construction = "closest";
% options.rom.bases_construction = "optimal";
% options.rom.bases_construction = "qr";

switch options.rom.bases_construction
    case "mthesis"
        X_hom = X_h - X_inhom;
        phi_hom = Om_POD(X_hom,M,options);
        phi_bc = POD(X_bc,Mbc);
        [phi_inhom,R_inhom,P] = get_phi_inhom(phi_bc,options);
    case "closest"
        [phi_h,weight_matrix] = Om_POD(X_h,M,options);
        phi_bc = get_velo_consis_phi_bc(X_bc,weight_matrix);
        [phi_inhom,R_inhom,P] = get_phi_inhom(phi_bc,options);
        phi_hom = homogeneous_projection(phi_h,options);
    case "optimal"
        [phi_h,weight_matrix] = Om_POD(X_h,M,options);
        phi_bc = get_velo_consis_phi_bc(X_bc,weight_matrix);
        [phi_inhom,R_inhom,P] = get_phi_inhom(phi_bc,options);
        X_hom = X_h - X_inhom;
        phi_hom = Om_POD(X_hom,M,options);
    case "qr"
        [phi_h,weight_matrix] = Om_POD(X_h,M,options);
        phi_bc = get_velo_consis_phi_bc(X_bc,weight_matrix);
        M_h = options.discretization.M;
        M_hphi = M_h*phi_h;
        rank_M_hphi = rank(M_hphi);
        [Q,R] = qr(M_hphi',0);
        Q1 = Q(:,1:rank_M_hphi);
        Q2 = Q(:,rank_M_hphi+1:end);
        R1 = R(1:rank_M_hphi,:);
        phi_hom = phi_h*Q2;
        phi_inhom = phi_h*Q1; % correct, but we also need R_inhom
%         [phi_inhom,R_inhom] = get_phi_inhom(phi_bc,options); % definitely wrong
        F_M = options.discretization.F_M;
        R_inhom = (R1*R1')\(R1*F_M*phi_bc);
        P = 1:size(phi_bc,2);
end
% testing
options.rom.bases_construction
figure; heatmap(phi_hom'*(Om.*phi_hom))
figure; heatmap(phi_inhom'*(Om.*phi_inhom))
norm(phi_hom'*(Om.*phi_inhom))
figure; heatmap(phi_hom'*(Om.*phi_inhom))
figure; heatmap(phi_bc'*phi_bc)
M_h = options.discretization.M;
F_M = options.discretization.F_M;
norm(M_h*phi_inhom*R_inhom-F_M*phi_bc(:,P))

B = phi_hom;
options.rom.B = B;
options.rom.phi_bc = phi_bc;
options.rom.phi_inhom = phi_inhom;
options.rom.R_inhom = R_inhom;

%% get Vbc into options (this has to be outside the j==1 if statement)
options.rom.Vbc = Vbc;

%% construct basis through SVD or eigenvalue problem
velocity_basis_construction;

%% check whether basis is divergence free
% this gives max div for each basis vector (columns of B):
% note that yM should not be included here, it was already used in
% subtracting Vbc from the snapshots matrix
div_basis = max(abs(options.discretization.M*B),[],1); %
% max over all snapshots:
maxdiv_basis = max(div_basis);
if (maxdiv_basis > 1e-14)
    warning(['ROM basis not divergence free: ' num2str(maxdiv_basis)]);
end

%% pressure recovery
pressure_basis_construction;

%% compute boundary condition approximation and inhomogeneous ROM basis
BC_and_lifting_function_basis_construction;

%% precompute ROM operators by calling operator_rom
% results are stored in options structure
if options.rom.bc_recon ~= 4 || options.rom.bc_recon ~= 2
    if (options.rom.precompute_convection == 1 || options.rom.precompute_diffusion == 1 || ...
            options.rom.precompute_force == 1 || options.rom.pressure_recovery == 1)
        disp('precomputing ROM operators...');
        precompute_start = toc;
        options = operator_rom(options);
        precompute_end(j) = toc-precompute_start
    end
end

%% initialize reduced order solution
ROM_initialization;

%%

if (options.rom.process_iteration_FOM == 1)
    % map back to velocity space to get statistics of initial velocity field
    % note that V will not be equal to the specified initial field, because
    % B*B' does not equal identity in general
    V  = getFOM_velocity(R,t,options);
    
    [maxdiv(1), umom(1), vmom(1), k(1)] = check_conservation(V,t,options);
    
    
    % overwrite the arrays with total solutions
    if (save_unsteady == 1)
        uh_total(n,:) = V(1:options.grid.Nu);
        vh_total(n,:) = V(options.grid.Nu+1:end);
    end
end

if (options.rom.pressure_recovery == 1)
    % get initial pressure that corresponds to the ROM velocity field
    q = pressure_additional_solve_ROM(R,t,options);
    p = getFOM_pressure(q,t,options);
    
    % overwrite the arrays with total solutions
    if (save_unsteady == 1)
        p_total(n,:)  = p;
    end
end


%% reduced order solution tests

% % test the offline implementation as follows:
% Rtest = rand(options.rom.M,1);
% Vtest = getFOM_velocity(Rtest,t,options);
%
% % with precomputing:
% conv_ROM_pre = convectionROM(Rtest,t,options,0);
% diff_ROM_pre = diffusionROM(Rtest,t,options,0);
%
% % without precomputing:
% [convu, convv, dconvu, dconvv] = convection(Vtest,Vtest,t,options,0);
% conv_ROM  = B'*(Om_inv.*[convu;convv]);
% [d2u,d2v,dDiffu,dDiffv] = diffusion(Vtest,t,options,0);
% diff_ROM  = B'*(Om_inv.*[d2u;d2v]);
%
% % compute error between the two versions
% error_conv_ROM = conv_ROM_pre - conv_ROM;
% error_diff_ROM = diff_ROM_pre - diff_ROM;
% plot(error_conv_ROM)
% hold on
% plot(error_diff_ROM);

%% load restart file if necessary
if (options.output.save_results==1)
    fprintf(fconv,'n            dt               t                res              maxdiv           umom             vmom             k\n');
    fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n',...
        n,dt,t,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
end

%% plot initial solution
if (rtp.show==1)
    run(rtp.file);
    % for movies, capture this frame
    if (rtp.movie==1)
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
end

%% start time stepping

dtn    = dt;

eps    = 1e-12;

disp('starting time-stepping...');

% addpath 'C:\Users\20201213\Documents\Uni\Master thesis\clean_code\INS2D\debug_stuff'
% addpath 'debug_stuff'
% jacobian_test_ROM

time_start = toc

coefficients = zeros(M,nt);

k_sum2 = zeros(nt,1);
k_delta = [];
k_analysis = [];

while(n<=nt)
    
    %% dynamic timestepping:
    % set_timestep;
    
    %%
    
    n = n+1;
    
    % time reversal (used in inviscid shear-layer testcase)
%     if (abs(t-t_end/2)<1e-12)
%         dt = -dt;
%     end
    
    % perform one time step with the time integration method
    if (method == 20)
            time_ERK_ROM;
    elseif (method == 21)
        if options.rom.bc_recon == 2
            pn = pressure_additional_solve(V,p,t,options);
            [R,p,iterations,k_delta,k_analysis] ...
                = time_IRK_ROM_notvelocityonly(R,p,t,dt,options,k_analysis);
%             time_IRK_ROM_notvelocityonly;
        else
            time_IRK_ROM;
        end
        
    else
        error('time integration method unknown');
    end
     
    coefficients(:,n) = R;
        
    % the velocities and pressure that are just computed are at
    % the new time level t+dt:
    t = t + dt;
    time(n) = t;
    
    % check residuals, conservation, write output files
    % this requires to go back to FOM level
    if (options.rom.process_iteration_FOM == 1)
        % map back from reduced space to full order model space
        % this is used for postprocessing purposes, e.g. evaluating the divergence
        % of the velocity field
        V = getFOM_velocity(R,t,options);

        if options.rom.bc_recon == 3
%             Np = options.grid.Np;
%             p = zeros(Np,1);
%             p = pressure_additional_solve(V,p,t,options); %wrong because
%             uses analytical derivative ydM
            yBCn = get_bc_vector_yBC(t-dt,options);
            yMn = get_yM(options,yBCn);
            yBCn1 = get_bc_vector_yBC(t,options);
            yMn1 = get_yM(options,yBCn1);
            p = pressure_additional_solve2(V,t,(yMn1-yMn)/dt,options);
            warning('only correct for forward Euler')
        elseif options.rom.bc_recon == 5
            Bp = options.rom.Bp;
            p = Bp*q;

            %% verbosity
%             yBCn = get_bc_vector_yBC(t-dt,options);
%             yMn = get_yM(options,yBCn);
%             yBCn1 = get_bc_vector_yBC(t,options);
%             yMn1 = get_yM(options,yBCn1);
%             p2 = pressure_additional_solve2(V,t,(yMn1-yMn)/dt,options);
%             norm(p-p2)
%             p_error = p-p2;
%             norm(B'*G*p_error)
            %%
        end
        
        process_iteration;
    end
    
    
end
disp('finished time-stepping...');
time_loop(j) = toc-time_start


V  = getFOM_velocity(R,t,options);