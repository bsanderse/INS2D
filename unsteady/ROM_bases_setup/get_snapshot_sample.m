function snapshot_sample = get_snapshot_sample(dt_sample,dt_snapshots,t_sample) 

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