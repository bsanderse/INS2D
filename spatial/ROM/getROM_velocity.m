function [R] = getROM_velocity(V,t,options)
% get ROM velocity coefficients

B   = options.rom.B;

if (options.rom.weighted_norm == 0)
    R   = B'*V;
elseif (options.rom.weighted_norm == 1)
    Om = options.grid.Om;
    R  = B'*(Om.*V);
end

