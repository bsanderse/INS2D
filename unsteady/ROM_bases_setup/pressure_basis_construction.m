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
%         cond_fac = 10^-6;
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
    
%     svd_end(j) = svd_end(j) + toc - svd_start2
    
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