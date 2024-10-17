function [V_svd,Vbc,snapshots] = load_snapshot_data(snapshot_data,options,dt_sample,t_sample)
        %% load snapshot data
        % assume that for parametric studies (e.g. changing number of modes, the
        % FOM data file does not change:
        % if (j==1)
        if true
            
            disp(['loading datafile...: ' snapshot_data]);
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
            
            if (snapshots.Re ~= options.fluid.Re)
                error('Reynolds numbers of snapshot data and current simulation do not match');
            end
            
            
            %% check whether snapshots are divergence free
            % this gives max div for each snapshot:
            div_snapshots = max(abs(options.discretization.M*V_total_snapshots + options.discretization.yM),[],1); %
            % max over all snapshots:
            maxdiv_snapshots = max(div_snapshots);
            if (maxdiv_snapshots > 1e-14)
                warning(['snapshots not divergence free: ' num2str(maxdiv_snapshots)]);
            end
            
            
            %% subtract non-homogeneous BC contribution:
            
            % note uh_total is stored as a Nt*Nu matrix, instead of the Nu*Nt matrix
            % which we use for the SVD
            Om     = options.grid.Om;
            Om_inv = options.grid.Om_inv;
            
            if (options.rom.rom_bc == 1)
                % check if the Vbc field has been stored as part of the FOM
                if (isfield(snapshots,'Vbc'))
                    Vbc = snapshots.Vbc;
                else
                    disp('computing Vbc field...');
                    f       = options.discretization.yM;
                    dp      = pressure_poisson(f,t,options);
                    Vbc     = - Om_inv.*(options.discretization.G*dp);
                end
                V_total_snapshots = V_total_snapshots - Vbc; % this velocity field satisfies M*V_total = 0
            else
                Vbc = zeros(Nu+Nv,1);
            end
            
            % sample dt can be used to get only a subset of the snapshots
            if (rem(dt_sample,snapshots.dt) == 0)
                % sample dt should be a multiple of snapshot dt:
                Nskip = dt_sample/snapshots.dt;
                % check if t_sample is multiple of dt_sample
                if (rem(t_sample,dt_sample) == 0)
                    Nsnapshots    = t_sample / snapshots.dt; %size(V_total,2);
                    snapshot_sample = 1:Nskip:Nsnapshots;
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