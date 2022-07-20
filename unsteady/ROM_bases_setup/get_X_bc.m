function X_bc = get_X_bc(t_start,t_end,dt,snapshot_sample,options)

ts_snapshots = t_start:dt:t_end;
ts_sample = ts_snapshots(snapshot_sample);

X_bc = zeros(length(get_bc_vector_yBC(0,options)),length(ts_sample));
for jj=1:length(ts_sample)
    t_j = ts_sample(jj);
    X_bc(:,jj) = get_bc_vector_yBC(t_j,options);
end