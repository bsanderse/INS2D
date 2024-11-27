if ~exist("IC_file",'var')
    IC_file = [options.case.project '_IC'];
end

if (exist(IC_file,'file'))
    
    t       = options.time.t_start;
    IC      = str2func(IC_file);    
    [u_start,v_start,p_start,options] = IC(t,options);
    
    V_start = [u_start(:);v_start(:)];

else
    
    error(['initial condition file ' IC_file ' not available']);
    
end
