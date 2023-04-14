t       = options.time.t_start;
if options.rom.initialize_with_snapshots
%     snapshot_data = 'RBC_data/matlab_data.mat';
    snapshot_data = '/export/scratch1/krishan/INS2D_scratch/data_generation/1e4/FOM/matlab_data.mat';
    disp(['loading datafile...: ' snapshot_data]);
    snapshots = load(snapshot_data,'uh_total','vh_total','p_total','T_total');
    
    u_start = snapshots.uh_total(1,:)';
    v_start = snapshots.vh_total(1,:)';
    p_start = snapshots.p_total(1,:)';
    T_start = snapshots.T_total(1,:)';
    V_start = [u_start(:);v_start(:)];
else
    
    file_name = [options.case.project '_IC'];
    
    
    if (exist(file_name,'file'))
        

        IC      = str2func(file_name);
        switch options.case.boussinesq
            case 'none'
                [u_start,v_start,p_start,options] = IC(t,options);
            case 'temp' %temperature equation added, need initial condition
                [u_start,v_start,p_start,T_start,options] = IC(t,options);
            otherwise
                error('wrong option provided for boussinesq');
        end
        V_start = [u_start(:);v_start(:)];
    else
        
        error(['initial condition file ' file_name ' not available']);
        
    end
    
end
%% This is done to compare the difference between best approximate
%   soluntion used as initial condition

%             snapshot_data = 'RBC_data/matlab_data.mat';
%             disp(['loading datafile...: ' snapshot_data]);
%             snapshots = load(snapshot_data,'u_best','v_best','p_total','T_best');
%             u_start = snapshots.u_best(options.grid.indu,1);
%             v_start = snapshots.v_best(options.grid.indv,1);
%             T_start = snapshots.T_best(:,2);
%             V_start = [u_start(:);v_start(:)];
