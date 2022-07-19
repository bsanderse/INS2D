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