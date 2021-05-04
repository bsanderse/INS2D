function [Vbc] = get_unsteadyVbc(t,options)
% get FOM pressure based on ROM coefficients

if options.rom.bc_recon == 1
    Vbc = options.rom.Vbc0*options.rom.abc(t);
else
    Om_inv = options.grid.Om_inv;
    
    options = set_bc_vectors(t,options);
    f       = options.discretization.yM;
    dp      = pressure_poisson(f,t,options);
    Vbc = - Om_inv.*(options.discretization.G*dp);
end