%  Main solver file for unsteady calculations with reduced order model


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
    
    % check input dimensions
    Nspace  = size(V_total_snapshots,1); % total number of unknowns (Nu+Nv) of the original model
    Nu      = options.grid.Nu;
    Nv      = options.grid.Nv;
    if (Nspace ~= Nu+Nv)
        error('The dimension of the snapshot matrix does not match the input dimensions in the parameter file');
    end
    
%     if (snapshots.Re ~= options.fluid.Re)
    if (snapshots.Re ~= Re)
        error('Reynolds numbers of snapshot data and current simulation do not match');
    end
    
    FOM_data_consistency_check(options,snapshots.options)
    
    
    %% check whether snapshots are divergence free
    if (options.BC.BC_unsteady==0)     % not working for unsteady BC
        % this gives max div for each snapshot:
        div_snapshots = max(abs(options.discretization.M*V_total_snapshots + options.discretization.yM),[],1); %
        % max over all snapshots:
        maxdiv_snapshots = max(div_snapshots);
        if (maxdiv_snapshots > 1e-14)
            warning(['snapshots not divergence free: ' num2str(maxdiv_snapshots)]);
        end
    end
    
    %% subtract non-homogeneous BC contribution:
    
    % note uh_total is stored as a Nt*Nu matrix, instead of the Nu*Nt matrix
    % which we use for the SVD
    Om     = options.grid.Om;
    Om_inv = options.grid.Om_inv;
    
    if (options.rom.bc_recon == 2) || (options.rom.bc_recon == 4) ...
            || (options.rom.bc_recon == 5  && options.verbosity.equivalence_cheat == 0)
        Vbc = zeros(Nu+Nv,1);
        snapshots.Vbc = Vbc;
    else %if (options.rom.bc_recon ~= 5)
        if (options.rom.rom_bc == 1)
            % check if the Vbc field has been stored as part of the FOM
            if (isfield(snapshots,'Vbc'))
                Vbc = snapshots.Vbc;
            else
                disp('computing Vbc field...');
                f       = options.discretization.yM;
                dp      = pressure_poisson(f,t,options);
                Vbc     = - Om_inv.*(options.discretization.G*dp);
                snapshots.Vbc = Vbc;
            end
            V_total_snapshots = V_total_snapshots - Vbc; % this velocity field satisfies M*V_total = 0
        elseif (options.rom.rom_bc == 2)
            if (isfield(snapshots,'Vbc'))
                Vbc = snapshots.Vbc;
            else
                %error('Vbc data not provided');
                disp('computing snapshot.Vbc')
%                 dt = snapshots.dt;
%                 t_end = snapshots.t_end;
                t_js = 0:dt:t_end;
                for jj=1:length(t_js)
                    t_j = t_js(jj);
                    options = set_bc_vectors(t_j,options);
                    f       = options.discretization.yM;
                    dp      = pressure_poisson(f,t_j,options);
                    Vbc(:,jj)  = - Om_inv.*(options.discretization.G*dp);
%                     Vbc(:,jj)  = Om_inv.*(options.discretization.G*dp); %pfusch
                    snapshots.Vbc = Vbc;
                end 
            end
            V_total_snapshots = V_total_snapshots - Vbc; % this velocity field satisfies M*V_total = 0
        else
            Vbc = zeros(Nu+Nv,1);
            snapshots.Vbc = Vbc;
        end
    end
    
    % sample dt can be used to get only a subset of the snapshots
    if (rem(dt_sample,dt_snapshots) == 0)
        % sample dt should be a multiple of snapshot dt:
        Nskip = dt_sample/dt_snapshots;
        % check if t_sample is multiple of dt_sample
        if (rem(t_sample,dt_sample) == 0)
            Nsnapshots    = round(t_sample / dt_snapshots); %size(V_total,2);
            snapshot_sample = 1:Nskip:(Nsnapshots+1);
        else
            error('sample dt is not an integer multiple of sample time');
        end
    else
        error('sample dt is not an integer multiple of snapshot dt');
    end
    
    % select snapshots
    V_svd = V_total_snapshots(:,snapshot_sample);
    
    clear V_total_snapshots;  
end

%% get Vbc into options (this has to be outside the j==1 if statement)
options.rom.Vbc = Vbc;


%% construct basis through SVD or eigenvalue problem
svd_start = toc;

% enforce momentum conservation (works for periodic domains)
if (options.rom.mom_cons == 1 && options.rom.weighted_norm == 0)
    
    e_u = zeros(Nspace,1);
    e_v = zeros(Nspace,1);
    e_u(1:Nu)     = 1;
    e_v(Nu+1:end) = 1;
    e = [e_u e_v];
    e = e / norm(e);
    
    % 1) construct (I-ee')*V_svd
    Vmod = V_svd - e*(e'*V_svd);
    % 2) take SVD
    [W,S,Z] = svd(Vmod,'econ');
    % 3) add e
    W = [e W];
    
    %     disp('error in representing vector y before truncating:');
    %     norm(Wext*Wext'*e - e,'inf')
    
elseif (options.rom.mom_cons == 1 && options.rom.weighted_norm == 1)
    
    Om_mat     = spdiags(Om,0,Nu+Nv,Nu+Nv);
    Om_sqrt    = spdiags(sqrt(Om),0,Nu+Nv,Nu+Nv);
    Om_invsqrt = spdiags(1./sqrt(Om),0,Nu+Nv,Nu+Nv);
    
    e_u = zeros(Nspace,1);
    e_v = zeros(Nspace,1);
    e_u(1:Nu)     = 1;
    e_v(Nu+1:end) = 1;
    e = [e_u e_v];
    % scale e such that e'*Om*e = I
    e = e / sqrt(norm(e'*(Om_mat*e)));
    
    % 1) construct (I-ee')*Om*V_svd
    Vmod = V_svd - e*(e'*(Om_mat*V_svd));
    % 2) apply weighting
    Vmod = Om_sqrt*Vmod;
    % 3) perform SVD
    [W,S,Z] = svd(Vmod,'econ');
    % 4) transform back
    W = Om_invsqrt*W;
    % 5) add e
    W = [e W];
    
elseif (options.rom.mom_cons == 0 && options.rom.weighted_norm == 0)
    
    % perform SVD
    %     [W,S,Z] = svd(V_svd,'econ');
    % getBasis can use different methods to get basis: SVD/direct/snapshot
    % method
    [W,S] = getBasis(V_svd,options);
    
elseif (options.rom.mom_cons == 0 && options.rom.weighted_norm == 1)
    
    Om_sqrt    = spdiags(sqrt(Om),0,Nu+Nv,Nu+Nv);
    Om_invsqrt = spdiags(1./sqrt(Om),0,Nu+Nv,Nu+Nv);
    
    % make weighted snapshot matrix
    Vmod = Om_sqrt*V_svd;
    % perform SVD
    %     [W,S,Z] = svd(Vmod,'econ');
    % getBasis can use different methods to get basis: SVD/direct/snapshot
    % method
    [W,S] = getBasis(Vmod,options);

    W0 = W;
    
    % transform back
    W = Om_invsqrt*W;

    %% mission: phi consistent phi bc
    norm(W0(:,1:M)-Vmod*Vmod'*W0(:,1:M)*diag(S(1:M).^-2))
    norm(W0(:,1:M)-Om_sqrt*V_svd*Vmod'*W0(:,1:M)*diag(S(1:M).^-2))
    norm(Om_invsqrt*W0(:,1:M)-V_svd*Vmod'*W0(:,1:M)*diag(S(1:M).^-2))

    if options.rom.rom_bc == 2
            t_js = t_start:dt:t_end;
    else
        t_js = 0;
        Mbc = 1;
    end
    X_bc = zeros(length(get_bc_vector_yBC(0,options)),length(t_js));
    for jj=1:length(t_js)
        t_j = t_js(jj);
        X_bc(:,jj) = get_bc_vector_yBC(t_j,options);
    end
    W_bc = X_bc*Vmod'*W0(:,1:Mbc)*diag(S(1:Mbc).^-2);

    phi_bc = W_bc(:,1:Mbc); % superfluous
    options.rom.phi_bc = phi_bc;
    for jj = 1:Mbc
        yBC = phi_bc(:,jj);
        Y_M(:,jj) = get_yM(options,yBC);
    end
    M_h = options.discretization.M;
    norm(M_h*W(:,1:Mbc)-Y_M) % not necesssarily true
    %%
    
else
    error('wrong option for weighted norm or momentum conservation');
    
end
% clear V_svd;

svd_end(j) = toc-svd_start

% take first M columns of W as a reduced basis
% maximum:
% M = size(Wu,2);
% reduction:
% M = floor(Nspace/100);
% M = 16;
M = options.rom.M;
% (better is to look at the decay of the singular values in S)
B  = W(:,1:M);
% Bu = Wu(:,1:M);
% Bv = Wv(:,1:M);
Bu = B(1:Nu,:);
Bv = B(Nu+1:end,:);
options.rom.B = B;
options.rom.Bu = Bu;
options.rom.Bv = Bv;
% options.rom.BuT = BuT;
% options.rom.BvT = BvT;
toc

% relative information content:
if (size(S,2)>1)
    Sigma = diag(S);
else
    Sigma = S;
end
% RIC  = sum(Sigma(1:M).^2)/sum(Sigma.^2);
% disp(['relative energy captured by SVD = ' num2str(RIC)]);
disp('commented some information out')

if (options.visualization.show_sigmas == 1)
    figure(123)
    semilogy(Sigma/Sigma(1),'s','displayname', 'singular values velocity snapshot matrix');
    hold on
end

% or alternatively
% semilogy(Sigma.^2/sum(Sigma.^2),'s');

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
if (options.rom.pressure_recovery == 1) || options.rom.bc_recon == 5
    disp('computing SVD of pressure snapshots...');
    svd_start2 = toc;
    % note p_total is stored as a Nt*Np matrix instead of Np*Nt which we use for
    % the SVD
    % use same snapshot_indx that was determined for velocity
    
    % select snapshots
    p_total_snapshots = snapshots.p_total';
    if (options.rom.pressure_mean == 1)    
        % subtract temporal mean
        options.rom.p_mean = mean(p_total_snapshots,2);
        p_total_snapshots = p_total_snapshots - options.rom.p_mean;
    end
    p_svd  = p_total_snapshots(:,snapshot_sample);

    % take first Mp columns of Wp as a reduced basis
    % (better is to look at the decay of the singular values in Sp to determine M)
    if (isfield(options.rom,'Mp'))
        Mp = options.rom.Mp;
    else
        % if not defined, use same number of modes as for velocity
        warning('number of pressure modes not defined, defaulting to number of velocity modes');
        Mp = options.rom.M;
    end

    if (options.rom.weighted_norm == 0)
        
        [Wp,Sp] = getBasis(p_svd,options,options.rom.Mp);
        
        % perform SVD
        %     [Wp,Sp,Zp] = svd(p_svd,'econ');
        
    elseif (options.rom.weighted_norm == 1)
        
        Np          = options.grid.Np;
        Omp         = options.grid.Omp;
        Omp_sqrt    = spdiags(sqrt(Omp),0,Np,Np);
        Omp_invsqrt = spdiags(1./sqrt(Omp),0,Np,Np);
        
        % make weighted snapshot matrix
        pmod = Omp_sqrt*p_svd;
        
        % getBasis can use different methods to get basis: SVD/direct/snapshot
        % method
        [Wp,Sp] = getBasis(pmod,options);
        
        % transform back
        Wp = Omp_invsqrt*Wp;
    end
    
    if options.rom.bc_recon == 5
        % for bc_recon == 5, existence of a solution to the ROM PPE requires in
        % general the invertibility of \hat L = \hat M \hat G <=>  \hat M has
        % full rank
        if M < Mp 
            warning('Sorry, pressure ROM basis is too big, is made smaller')
            Mp = M;
        end
        M_h = options.discretization.M;
        Bp = Wp(:,1:Mp);
%         while rank(Bp'*M_h*B)<Mp % prone to machine precision problems
%         rank = sum(abs(svd(Bp'*M_h*B))>10^-10);
%         cond_fac = 10^-8;
        cond_fac = 10^-6;
        sing_vals = svd(Bp'*M_h*B);
        rank_ = sum(abs(sing_vals)>cond_fac*max(abs(sing_vals))); % avoid badly scaled hatL
        while rank_<Mp
            warning('Sorry, pressure ROM basis is too big, is made smaller')
            Mp = rank_;
            Bp = Wp(:,1:Mp);
%             rank = sum(abs(svd(Bp'*M_h*B))>10^-10);
            sing_vals = svd(Bp'*M_h*B);
            rank_ = sum(abs(sing_vals)>cond_fac*max(abs(sing_vals))); % avoid badly scaled hatL
        end
    else
        Bp = Wp(:,1:Mp);
    end 
    Mp
    options.rom.Bp = Bp;
    
    svd_end(j) = svd_end(j) + toc - svd_start2
    
    hold on
    if (size(Sp,2)>1)
        SigmaP = diag(Sp);
    else
        SigmaP = Sp;
    end
    if (options.visualization.show_sigmas == 1)
        semilogy(SigmaP/SigmaP(1),'o','displayname', 'singular values pressure snapshot matrix');
    end
end
if (options.visualization.show_sigmas == 1)
    ylabel("\sigma_i/\sigma_1")
    xlabel("mode index")
    title('singular values')
    legend('show')
    if (exist('fig_destination') && j==Nsim)
        savefig([fig_destination '/singular values'])
    end
end

%% compute boundary condition approximation and inhomogeneous ROM basis
if (options.rom.bc_recon == 3) || (options.rom.bc_recon == 5) 
%     if options.rom.rom_bc == 2
%         dt = snapshots.dt;
%         t_end = snapshots.t_end;

%         if options.rom.rom_bc == 2
%             t_js = t_start:dt:t_end;
%         else
%             t_js = 0;
%             Mbc = 1;
%         end
%         X_bc = zeros(length(get_bc_vector_yBC(0,options)),length(t_js));
%         for jj=1:length(t_js)
%             t_j = t_js(jj);
%             X_bc(:,jj) = get_bc_vector_yBC(t_j,options);
        end
        if options.rom.rom_bc == 2
            [U_bc,S_bc,V_bc] = svd(X_bc,'econ');
            if (options.visualization.show_sigmas == 1)
                Sigma_bc = diag(S_bc);
                figure
                semilogy(Sigma_bc/Sigma_bc(1),'s','displayname', 'singular values Vbc snapshot matrix');
                ylabel("\sigma_i/\sigma_1")
                xlabel("mode index")
                title('singular values')
                legend('show')
                if (exist('fig_destination') && j==Nsim)
                    savefig([fig_destination '/singular values'])
                end
            end
        else
            U_bc = X_bc/norm(X_bc);
        end
        
%         cond_fac = 10^-6;
%         X_bc_rank = sum(abs(Sigma_bc/Sigma_bc(1))>cond_fac);
%         if Mbc > X_bc_rank
%             Mbc = X_bc_rank;
%         end
% 
%         phi_bc = U_bc(:,1:Mbc);
%         options.rom.phi_bc = phi_bc;
%         if (options.rom.bc_recon == 3) || options.verbosity.equivalence_cheat == 1
%             for jj = 1:Mbc
%                 yBC = phi_bc(:,jj);
%                 Y_M(:,jj) = get_yM(options,yBC);
%             end
            L = options.discretization.A;
            Gx   = options.discretization.Gx;
            Gy   = options.discretization.Gy;
            G = [Gx;Gy];
            Om = options.grid.Om;
            Om_inv = options.grid.Om_inv;
            %     tilde_phi_inhom = Om_inv.*(G*(L\Y_M));
            tilde_phi_inhom = -Om_inv.*(G*(L\Y_M)); %pfusch

            [Q_inhom,R_inhom] = qr(sqrt(Om).*tilde_phi_inhom,0); %alternative: take first vec of tilde phi inhom
            M_inhom = rank(tilde_phi_inhom);
            Q_1_inhom = -Q_inhom(:,1:M_inhom);
            R_inhom = -R_inhom(1:M_inhom,:);
            phi_inhom = sqrt(Om_inv).*Q_1_inhom;

            options.rom.phi_inhom = phi_inhom;
            options.rom.R_inhom = R_inhom;
%%
%             warning('phi_inhom computation manipulated')
%             M_h = options.discretization.M;
%             [Qq,Rr] = qr((M_h*phi_inhom)');
%             rankk = rank(M_h*phi_inhom);
%             Qq1 = Qq(:,1:rankk);
%             phi_inhom_star = phi_inhom*Qq1;
%             R_inhom_star = Qq1'*R_inhom;
% 
%             options.rom.phi_inhom = phi_inhom_star;
%             options.rom.R_inhom = R_inhom_star;
%%
        else
            options.rom.phi_inhom = 0;
            options.rom.R_inhom = 0;
        end
end

if options.verbosity.equivalence_cheat == 1
    B = [B phi_inhom];
    options.rom.B = B;
    M = M + M_inhom;
    options.rom.M = M;


    %%
if options.rom.bc_recon == 5
        % for bc_recon == 5, existence of a solution to the ROM PPE requires in
        % general the invertibility of \hat L = \hat M \hat G <=>  \hat M has
        % full rank
        if M < Mp 
            warning('Sorry, pressure ROM basis is too big, is made smaller')
            Mp = M;
        end
        M_h = options.discretization.M;
        Bp = Wp(:,1:Mp);
%         while rank(Bp'*M_h*B)<Mp % prone to machine precision problems
%         rank = sum(abs(svd(Bp'*M_h*B))>10^-10);
%         cond_fac = 10^-8;
        cond_fac = 10^-6;
        sing_vals = svd(Bp'*M_h*B);
        rank_ = sum(abs(sing_vals)>cond_fac*max(abs(sing_vals))); % avoid badly scaled hatL
        while rank_<Mp
            warning('Sorry, pressure ROM basis is too big, is made smaller')
            Mp = rank_;
            Bp = Wp(:,1:Mp);
%             rank = sum(abs(svd(Bp'*M_h*B))>10^-10);
            sing_vals = svd(Bp'*M_h*B);
            rank_ = sum(abs(sing_vals)>cond_fac*max(abs(sing_vals))); % avoid badly scaled hatL
        end
    else
        Bp = Wp(:,1:Mp);
    end 
    options.rom.Bp = Bp;
    %%
end

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
% we expand the part of the solution vector that is div-free in terms of
% B*R
% V = B*R + Vbc

if options.rom.bc_recon ~= 5
% get the coefficients of the ROM
R = getROM_velocity(V,t,options);

% for projected-divergence-free ROM, enforce projected divergence-freeness
% if options.rom.bc_recon == 5
else
    % construct ROM divergence operator
        Bp = options.rom.Bp;
        hatM = Bp'*options.discretization.M*B;
        options.rom.hatM = hatM;

%     Bp = options.rom.Bp;
%     hatM = options.rom.hatM;
    hatG = -hatM';
    hatL = hatM*hatG;
    yM = -options.discretization.yM;
    hatyM = Bp'*yM;
    %     phi_bc = options.rom.phi_bc;
    %     yM = phi_bc*phi_bc'*yM;

    %     bstar = hatL\(hatM*R-Bp'*yM);
    %     Rstar = R-hatG*bstar;
    %     R = Rstar;

    % QR based method
    Om = options.grid.Om;
%     [Q_,R_] = qr(hatM');
%     Q_1 = Q_(:,1:Mp);
%     Q_2 = Q_(:,Mp+1:end);
%     R_1 = R_(1:Mp,1:Mp);
%     R = Q_1*hatyM/(R_1') + Q_2*Q_2'*B'*(Om.*V);
% 
%         %test
%     norm(hatM*R-hatyM)

    % second try
    M_h = options.discretization.M;
    [Q__,R__] = qr((M_h*B)');
    H = rank(M_h*B);
    Q_1_ = Q__(:,1:H);
    Q_2_ = Q__(:,H+1:end);
    R_1_ = R__(1:H,:);
    a_1 = (R_1_*R_1_')\(R_1_*yM);
    a_2 = Q_2_'*B'*(Om.*V);
    R = Q_1_*a_1 + Q_2_*a_2;

    % tests
    norm(M_h*B*Q_1_*a_1-yM)
    diff_ = B*R-V;
    norm(diff_)
    norm(M_h*diff_)
    norm(a_2)
    norm(diff_'*(Om.*V))
    diff_2 = B*Q_1_*a_1 - V;

%     [QQ,RR] = qr(full(M_h)); % expensive!
%     HH = rank(M_h);
%     QQ1 = QQ(:,1:HH);

%     Q1t = options.discretization.Q1t;
%     Q2t = options.discretization.Q2t;
%     R1 = options.discretization.R1;
% 
%     a1 = (R1*R1')\(R1*yM); % only valid if ROM basis includes all relevant inhomogeneous modes
%     a11 = Q1t'*(Om.*V); % only valid if ROM basis includes all relevant inhomogeneous modes
% 
%     norm(a1-a11)
% 
%     a0 = B'*(Om.*V);
%     norm(M_h*B*a0)
%     
%     Qt = [Q1t Q2t];
%     T = Qt'*(Om.*B);
%     a2 = Q2t'*(Om.*V);
%     at = [a1; a2];
%     a = (T'*T)\(T'*at); % does not yield exact solution to T*a = at, but approximation
%     norm(T*a-at)
%     norm(R-a)
%     R = a;

%     Vr = B*a;
%     norm(M_h*Vr-yM)
%     norm(Vr-V)
end

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