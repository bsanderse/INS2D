function [R] = getROM_velocity2(V,t,options,basis)
% get ROM velocity coefficients

% B   = options.rom.B;
B = basis;

% subtract boundary condition contribution (zero if not used)
% if V is a NV*Nt matrix, then this vector is subtracted from each column
if options.rom.bc_recon ~= 5

    if options.rom.bc_recon
        if (options.rom.rom_bc == 2)
            V = V - get_unsteadyVbc(t,options); %not sure if necessary as Vbc should be orthogonal to B
        else
            V   = V - options.rom.Vbc(:,1);
        end
    end

end

if (options.rom.weighted_norm == 0)
    R   = B'*V;
elseif (options.rom.weighted_norm == 1)
    Om = options.grid.Om;
    R  = B'*(Om.*V);
end
