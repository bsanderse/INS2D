file_name = [options.case.project '_IC'];

if (exist(file_name,'file'))
    
    t       = options.time.t_start;
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

%             snapshot_data = '/export/scratch1/krishan/INS2D_scratch/perturbation_study/ROM_data_used_for_initial_condition_forFOM/matlab_data.mat';
%             disp(['loading datafile...: ' snapshot_data]);
%             snapshots = load(snapshot_data,'V_best','T_best');
%             u_start = snapshots.V_best(options.grid.indu,1);
%             v_start = snapshots.V_best(options.grid.indv,1);
%             T_start = snapshots.T_best(:,2);
%             V_start = [u_start(:);v_start(:)];
            
% /export/scratch1/krishan/INS2D_scratch/perturbation_study/ROM_bestapprox_initial_solution