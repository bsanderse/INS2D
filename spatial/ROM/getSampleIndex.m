function index = getSampleIndex(rom)

% sample dt can be used to get only a subset of the snapshots
dt_sample    = rom.dt_sample;
dt_snapshots = rom.dt_snapshots;
t_sample     = rom.t_sample;
if (rem(dt_sample,dt_snapshots) == 0)
    % sample dt should be a multiple of snapshot dt:
    Nskip = dt_sample/dt_snapshots;
    % check if t_sample is multiple of dt_sample
    if (rem(t_sample,dt_sample) == 0)
        Nsnapshots    = t_sample / dt_snapshots; %size(V_total,2);
        index = 1:Nskip:Nsnapshots;
    else
        error('sample dt is not an integer multiple of sample time');
    end
else
    error('sample dt is not an integer multiple of snapshot dt');
end