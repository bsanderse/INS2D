%  Main solver file for unsteady calculations with reduced order model

addpath('unsteady/ROM_bases_setup/');

%% load snapshot data
if (j==1) || changing_snapshotdata
    disp('loading data...');
    snapshots = load(snapshot_data,'uh_total','vh_total','p_total','dt','t_end','Re','k','umom','vmom','maxdiv','Vbc','options' ...
       ... ,'k_diff', 'k_conv', 'k_pres', 'k_force', 'k_diffBC', 'k_presBC','k_obc'
        );
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
    
    %% analysis
%     F_M = options.discretization.F_M;
%     [U_M,S_M,V_M] = svd(F_M*X_bc);
%     Sigma_M = diag(S_M);
%     figure
%     semilogy(Sigma_M/Sigma_M(1),'s','displayname', 'singular values of y_M snapshots');
%     legend('show')
    %%
    
    %% analysis
%     [U_h,S_h,V_h] = svd(X_h);
%     Sigma_h = diag(S_h);
%     figure
%     semilogy(Sigma_h/Sigma_h(1),'s','displayname', 'singular values of V_h snapshots');
%     legend('show')
    %%
end

%% get Vbc into options (this has to be outside the j==1 if statement)
options.rom.Vbc = Vbc;

%%

% cond_fac = 1e-10;
% cond_fac = 1e-20;
cond_fac = 1e-20;


if true
% if false
    
% options.rom.bases_construction = "mthesis";
% options.rom.bases_construction = "closest";
% % options.rom.bases_construction = "optimal";
% options.rom.bases_construction = "qr";

options.rom.Mbc = Mbc; % botch


switch options.rom.bases_construction
    case "mthesis"
        X_hom = X_h - X_inhom;
        [phi_hom,~,M] = Om_POD(X_hom,M,options,cond_fac);
        [phi_bc,Mbc] = POD(X_bc,Mbc,cond_fac);
        [phi_inhom,R_inhom,P,tilde_phi_inhom1] = get_phi_inhom(phi_bc,options);
    case "closest"
        [phi_h,weight_matrix,M] = Om_POD(X_h,M,options,cond_fac);
        phi_bc = get_velo_consis_phi_bc(X_bc,weight_matrix,false);
        [phi_inhom,R_inhom,P] = get_phi_inhom(phi_bc,options);
        phi_hom = homogeneous_projection(phi_h,options);
    case "optimal"
        [phi_h,weight_matrix,M] = Om_POD(X_h,M,options,cond_fac);
        phi_bc = get_velo_consis_phi_bc(X_bc,weight_matrix,false);
        [phi_inhom,R_inhom,P] = get_phi_inhom(phi_bc,options);
        X_hom = X_h - X_inhom;
        [phi_hom,~,M] = Om_POD(X_hom,M,options,cond_fac);
    case "qr"
        [phi_h,weight_matrix,M] = Om_POD(X_h,M,options,cond_fac);
        phi_bc = get_velo_consis_phi_bc(X_bc,weight_matrix,false);
        M_h = options.discretization.M;
        M_hphi = M_h*phi_h;
        rank_M_hphi = rank(M_hphi);
        [Q,R] = qr(M_hphi',0);
        %% experiment -> didn't help
%         [Q,R,P] = qr(M_hphi',0);
%         % test
%         norm((Q*R)'-M_hphi(P,:))
        %%
        Q1 = Q(:,1:rank_M_hphi);
        Q2 = Q(:,rank_M_hphi+1:end);
        R1 = R(1:rank_M_hphi,:);
        %% next experiment -> didn't help; presumption: would require to update phi_bc as well, don't know how that could be done
%         cond_fac = 1e-6;
%         relevance = (sum(abs(Q1'*M_hphi'),2)>cond_fac);
%         indices = find(relevance);
%         Q1 = Q1(:,indices);
%         R1 = R(indices,:);
        %%
        phi_hom = phi_h*Q2;
        phi_inhom = phi_h*Q1; % correct, but we also need R_inhom
%         [phi_inhom,R_inhom] = get_phi_inhom(phi_bc,options); % definitely wrong
        F_M = options.discretization.F_M;
        R_inhom = (R1*R1')\(R1*F_M*phi_bc);
        R_inhom = - R_inhom; % botch
        P = 1:size(phi_bc,2);
        
        condition1 = cond(R1*R1')
end
%testing
% options.rom.bases_construction
% figure; heatmap(phi_hom'*(Om.*phi_hom))
% figure; heatmap(phi_inhom'*(Om.*phi_inhom))
% norm(phi_hom'*(Om.*phi_inhom))
% figure; heatmap(phi_hom'*(Om.*phi_inhom))
% figure; heatmap(phi_bc'*phi_bc)
% M_h = options.discretization.M;
% F_M = options.discretization.F_M;
% % norm(M_h*phi_inhom*R_inhom-F_M*phi_bc(:,P))
conditions2(j) = cond(R_inhom)

if options.rom.bc_recon == 3

B = phi_hom;

% options.rom.B1 = B;
% options.rom.phi_bc1 = phi_bc;
% options.rom.phi_inhom1 = phi_inhom;
% options.rom.R_inhom1 = R_inhom;

options.rom.phi_inhom = phi_inhom;
options.rom.R_inhom = R_inhom;
M_inhom = size(phi_inhom,2);
options.rom.M_inhom = M_inhom;

elseif options.rom.bc_recon == 5
    switch options.rom.bases_construction
        case "mthesis"
            B = [phi_hom phi_inhom];
        case "closest"
            B = [phi_hom phi_inhom];
            
% %             F_M = options.discretization.F_M;
% %             Bp = orthonormalize(F_M*phi_bc); 
%             % in theory, F_M*phi_bc = M_h*phi_h
%             % in practice, however, rank computations are error-prone
%             M_h = options.discretization.M;
%             Bp = orthonormalize(M_h*phi_inhom);
%             % M_h*phi_inhom should span the same space as F_M*phi_bc while
%             % having smaller or equal (computed) rank
        case "optimal"
            B = [phi_hom phi_inhom];
        case "qr"
            B = phi_h;
%             Bp = orthonormalize(R1');
    end
    
%     options.rom.Bp = Bp;
%     Mp = size(Bp,2);
%     options.rom.Mp = Mp;

condition6(j) = cond(options.discretization.M*B)
condition7(j) = cond(options.discretization.M*phi_inhom)

end

options.rom.B = B;
M = size(B,2);
options.rom.M = M;

options.rom.phi_bc = phi_bc;
Mbc = size(phi_bc,2);
options.rom.Mbc = Mbc;
end


%%

% M = options.rom.M; %botch
% Mbc = options.rom.Mbc; %botch


% %% construct basis through SVD or eigenvalue problem
% velocity_basis_construction;
% 
% %% check whether basis is divergence free
% % this gives max div for each basis vector (columns of B):
% % note that yM should not be included here, it was already used in
% % subtracting Vbc from the snapshots matrix
% div_basis = max(abs(options.discretization.M*B),[],1); %
% % max over all snapshots:
% maxdiv_basis = max(div_basis);
% if (maxdiv_basis > 1e-14)
%     warning(['ROM basis not divergence free: ' num2str(maxdiv_basis)]);
% end

%% pressure recovery
if options.rom.bases_construction == "qr"
%     M_h = options.discretization.M;
%     Bp = orthonormalize(M_h*phi_inhom,false);
%     options.rom.Bp = Bp;
%     
%     Bp = orthonormalize(R1');
%     options.rom.Bp2 = Bp;
    
    Bp = orthonormalize(R1'/(R1*R1'));
%     options.rom.Bp3 = Bp;
    
    options.rom.Bp = Bp;
else
%     pressure_basis_construction;

    M_h = options.discretization.M;
    Bp = orthonormalize(M_h*phi_inhom,false);
    options.rom.Bp = Bp;
end

%% compute boundary condition approximation and inhomogeneous ROM basis
% BC_and_lifting_function_basis_construction;


%% testing/comparing
% B = options.rom.B;
% B1 = options.rom.B1;
% norm(B-B1)
% phi_inhom = options.rom.phi_inhom;
% phi_inhom1 = options.rom.phi_inhom1;
% norm(phi_inhom-phi_inhom1)
% norm(phi_inhom+phi_inhom1)
% phi_bc = options.rom.phi_bc;
% phi_bc1 = options.rom.phi_bc1;
% norm(phi_bc-phi_bc1)
% R_inhom = options.rom.R_inhom;
% R_inhom1 = options.rom.R_inhom1;
% norm(R_inhom-R_inhom1)
% norm(R_inhom+R_inhom1)
% 
% B = B1;
% options.rom.B = B;
% phi_bc = phi_bc1;
% options.rom.phi_bc = phi_bc;

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
%             warning('only correct for forward Euler')
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