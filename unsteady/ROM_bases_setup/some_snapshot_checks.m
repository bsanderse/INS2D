%%

    if (Nspace ~= Nu+Nv)
        error('The dimension of the snapshot matrix does not match the input dimensions in the parameter file');
    end
    
%     if (snapshots.Re ~= options.fluid.Re)
    if (snapshots.Re ~= Re)
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