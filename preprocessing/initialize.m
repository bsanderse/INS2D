file_name = [options.case.project '_IC'];

if (exist(file_name,'file'))
    
    t       = options.time.t_start;
    IC      = str2func(file_name);    
    [u_start,v_start,p_start,options] = IC(t,options);
    
else
    
    error(['initial condition file ' file_name ' not available']);
    
end
