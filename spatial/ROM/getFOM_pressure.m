function [p] = getFOM_pressure(q,t,options)
% get FOM pressure based on ROM coefficients

Bp  = options.rom.Bp;
p   = Bp*q;
if (options.rom.pressure_mean == 1)
    p = p + options.rom.p_mean;
end

