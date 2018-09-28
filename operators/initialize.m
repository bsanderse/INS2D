file_name = [options.case.project '_IC'];

if (exist(file_name,'file'))
    
    IC    = str2func(file_name);    
    [u,v] = IC(t,options);

else
    
    error(['initial condition file ' file_name ' not available']);
    
end
