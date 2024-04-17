function V_snapshots = load_snapshot_data(snapshot_data,dt_sample,t_sample)

disp('loading data...');
    snapshots = load(snapshot_data,'uh_total','vh_total','p_total','dt','t_end','Re','k','umom','vmom','maxdiv','Vbc','options' );
    
    % dt that was used for creating the snapshot matrix:
    dt_snapshots = snapshots.dt;
    % velocity field as snapshot matrix:
    V_total_snapshots = [snapshots.uh_total';snapshots.vh_total'];
     
    snapshot_sample = get_snapshot_sample(dt_sample,dt_snapshots,t_sample);
    
    % select snapshots
    V_snapshots = V_total_snapshots(:,snapshot_sample);
