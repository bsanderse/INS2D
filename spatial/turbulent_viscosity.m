function  nu_t = turbulent_viscosity(S_abs,options)

visc = options.case.visc;

% compute turbulent viscosity based on S_abs
switch visc
    case 'LES' % Smagorinsky
        
        C_S  = options.visc.Cs;
        filter_length = deltax; % =sqrt(FV size) for uniform grids
        
        nu_t = (C_S^2)*(filter_length^2)*S_abs;
        
    case 'qr'  % q-r
        % q-r
        C_d  = deltax^2/8;
        nu_t = C_d * 0.5 * S_abs * (1 - alfa / C_d)^2;
        
    case 'ML' % mixing-length
        lm   = options.visc.lm; % mixing length
        nu_t = (lm^2)*S_abs;
        
    otherwise
        error('wrong value for visc parameter');
        
end