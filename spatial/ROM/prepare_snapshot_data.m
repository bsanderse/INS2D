function [V_svd,Vbc,snapshots,Nspace] = prepare_snapshot_data(snapshot_data,options,multiple_train_data)

if multiple_train_data == 1
    V_svd = [];
    Vbc = [];
    for j = numel(snapshot_data)
    [V_svd_j,Vbc_j,snapshots_j,Nspace] = prepare_snapshot_data(char(snapshot_data(j)),options,0);
    V_svd = [V_svd, V_svd_j];
    Vbc = [Vbc, Vbc_j];
    end
    snapshots = snapshots_j;
else
    disp('loading data...');
    snapshots = load(snapshot_data,'uh_total','vh_total','p_total','dt','t_end','Re','k','umom','vmom','maxdiv','Vbc');
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
    
    if (snapshots.Re ~= options.fluid.Re && options.rom.vary_re==0)
        error('Reynolds numbers of snapshot data and current simulation do not match');
    end
    
    
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
    
    if (options.rom.bc_recon == 2)
        Vbc = 0*Om;
        snapshots.Vbc = Vbc;
    else
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
                error('Vbc data not provided');
            end
            V_total_snapshots = V_total_snapshots - Vbc; % this velocity field satisfies M*V_total = 0
        else
            Vbc = zeros(Nu+Nv,1);
            snapshots.Vbc = Vbc;
        end
    end
    
    dt_sample = options.rom.dt_sample;
    t_sample  = options.rom.t_sample;
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
