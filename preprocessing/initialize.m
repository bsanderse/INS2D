if file_format == 1
    file_name = 'IC_';
else
    file_name = [options.case.project '_IC'];
end

if (exist(file_name,'file'))
    
    t       = options.time.t_start;
    IC      = str2func(file_name);    
    [u_start,v_start,p_start,options] = IC(t,options);
    V_start = [u_start(:);v_start(:)];
else
    
    error(['initial condition file ' file_name ' not available']);
    
end
