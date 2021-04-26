function [R] = getROM_velocity(V,t,options)
% get ROM velocity coefficients

B   = options.rom.B;

% subtract boundary condition contribution (zero if not used)
% if V is a NV*Nt matrix, then this vector is subtracted from each column
if (options.rom.rom_bc == 2)
    V = V - get_unsteadyVbc(t,options);
else
    V   = V - options.rom.Vbc;
end

if (options.rom.weighted_norm == 0)
    R   = B'*V;
elseif (options.rom.weighted_norm == 1)
    Om = options.grid.Om;
    R  = B'*(Om.*V);
end
