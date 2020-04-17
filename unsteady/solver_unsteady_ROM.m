%  Main solver file for unsteady calculations with reduced order model


%% load snapshot data

disp('loading data and making SVD...');
snapshots = load(snapshot_data,'uh_total','vh_total','p_total','dt','t_end','Re','k','umom','vmom','maxdiv');
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

if (snapshots.Re ~= Re)
    error('Reynolds numbers of snapshot data and current simulation do not match');
end


%% construct economic SVD
% note uh_total is stored as a Nt*Nu matrix, instead of the Nu*Nt matrix 
% which we use for the SVD
Om     = options.grid.Om;
Om_inv = options.grid.Om_inv;

% subtract non-homogeneous BC contribution
if (options.rom.rom_bc == 1)
    f       = options.discretization.yM;
    dp      = pressure_poisson(f,t,options);
    Vbc     = - Om_inv.*(options.discretization.G*dp);
    V_total_snapshots = V_total_snapshots - Vbc; % this velocity field satisfies M*V_total = 0
else
    Vbc = zeros(Nu+Nv,1);
end
options.rom.Vbc = Vbc;



% sample dt is multiple of snapshot dt:
if (rem(dt_sample,dt_snapshots) == 0)
    Nskip = dt_sample/dt_snapshots;
    % check if t_sample is multiple of dt_sample
    if (rem(t_sample,dt_sample) == 0)
        Nsnapshots    = t_sample / dt_snapshots; %size(V_total,2);
        snapshot_indx = 1:Nskip:Nsnapshots;
    else
        error('sample dt is not an integer multiple of sample time');
    end
else
    error('sample dt is not an integer multiple of snapshot dt');
end


% select snapshots
V_svd = V_total_snapshots(:,snapshot_indx);

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
    [W,S,Z] = svd(V_svd,'econ');
    
elseif (options.rom.mom_cons == 0 && options.rom.weighted_norm == 1)
    
    Om_sqrt    = spdiags(sqrt(Om),0,Nu+Nv,Nu+Nv);
    Om_invsqrt = spdiags(1./sqrt(Om),0,Nu+Nv,Nu+Nv);
    
    % make weighted snapshot matrix
    Vmod = Om_sqrt*V_svd;
    % perform SVD
    [W,S,Z] = svd(Vmod,'econ');
    % transform back
    W = Om_invsqrt*W;
    
else
    error('wrong option for weighted norm or momentum conservation');
    
end

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
Sigma = diag(S);
RIC  = sum(Sigma(1:M).^2)/sum(Sigma.^2);
disp(['relative energy captured by SVD = ' num2str(RIC)]);
figure
semilogy(Sigma/Sigma(1),'s');
% or alternatively
% semilogy(Sigma.^2/sum(Sigma.^2),'s');

%% pressure recovery
if (options.rom.pressure_recovery == 1)
    disp('computing SVD of pressure snapshots...');

    % note p_total is stored as a Nt*Np matrix instead of Np*Nt which we use for
    % the SVD
    % use same snapshot_indx that was determined for velocity
    
    % select snapshots
    p_total_snapshots = snapshots.p_total';
    p_svd  = p_total_snapshots(:,snapshot_indx);
    
    % perform SVD
    [Wp,Sp,Zp] = svd(p_svd,'econ');
    
    % take first Mp columns of Wp as a reduced basis
    % (better is to look at the decay of the singular values in Sp to determine M)
    if (isfield(options.rom,'Mp'))
        Mp = options.rom.Mp;
    else
        % if not defined, use same number of modes as for velocity
        warning('number of pressure modes not defined, defaulting to number of velocity modes');
        Mp = options.rom.M;
    end
    Bp = Wp(:,1:Mp);
    options.rom.Bp = Bp;
    
    toc
    
    hold on
    SigmaP = diag(Sp);
    semilogy(SigmaP/SigmaP(1),'o');
       
end

%% precompute ROM operators by calling operator_rom
% results are stored in options structure
disp('precomputing ROM operators...');
options = operator_rom(options);


%% initialize reduced order solution
% we expand the part of the solution vector that is div-free in terms of
% B*R
V  = V - Vbc; % subtract boundary condition contribution (zero if not used)

% get the coefficients of the ROM
R = getROM_velocity(V,t,options);

% map back to velocity space to get statistics of initial velocity field
% note that V will not be equal to the specified initial field, because
% B*B' does not equal identity in general
V  = getFOM_velocity(R,t,options);

[maxdiv(1), umom(1), vmom(1), k(1)] = check_conservation(V,t,options);

if (options.rom.pressure_recovery == 1)
    % get initial pressure that corresponds to the ROM velocity field
    q = pressure_additional_solve_ROM(R,t,options);
    p = getFOM_pressure(q,t,options);
end

% overwrite the arrays with total solutions
if (steady==0 && save_unsteady == 1)
    uh_total(n,:) = V(1:options.grid.Nu);
    vh_total(n,:) = V(options.grid.Nu+1:end);
    p_total(n,:)  = p;
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

time_start = toc

% while (abs(t)<=(t_end-dt+eps))
% rev = 0;
while(n<=nt)
    
    %% dynamic timestepping:
    % set_timestep;
    
    %%
    
    n = n+1;
       
    % perform one time step with the time integration method
    if (method == 20)
        time_ERK_ROM;
    elseif (method == 21)
        time_IRK_ROM;
    else
        error('time integration method unknown');
    end
    
    
    % the velocities and pressure that are just computed are at
    % the new time level t+dt:
    t = t + dt;
    time(n) = t;
    
    % check residuals, conservation, write output files
    % this requires to go back to FOM level
    if (options.rom.process_iteration_FOM == 1)
        process_iteration;
    end
    
    
end
disp('finished time-stepping...');
time_loop = toc-time_start