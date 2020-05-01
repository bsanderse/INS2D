function [q] = getROM_pressure(p,t,options)
% get ROM pressure coefficients

Bp   = options.rom.Bp;

if (options.rom.weighted_norm == 0)
    q   = Bp'*p;
elseif (options.rom.weighted_norm == 1)
    Omp = options.grid.Omp;
    q   = Bp'*(Omp.*p);
end

