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
