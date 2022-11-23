function [R,q] = initializeROM(V,p,t,options)
% compute ROM variables from FOM variables

% we expand the part of the solution vector that is div-free in terms of
% B*R
% V = B*R + Vbc
% V_orig = V;
% get the coefficients of the ROM
R = getROM_velocity(V,t,options);


if (options.rom.div_free == 0)
    % get initial ROM pressure by projecting initial presssure
    q = getROM_pressure(p,t,options);
    
    % make initial condition 'divergence free'
    f  = options.rom.Mdiv*R + options.rom.yMt(:,1);
    dq = pressure_poisson_ROM(f,t,options);
    R  = R - options.rom.G*dq;
    
elseif (options.rom.div_free == 1)
    if (options.rom.pressure_recovery == 1)
        % get initial pressure that corresponds to the ROM velocity field
        q = pressure_additional_solve_ROM(R,t,options);        
    else
        % q is not used, set to arbitrary value
        q = 0;
    end
end