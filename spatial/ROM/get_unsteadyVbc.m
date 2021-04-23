function [Vbc] = get_unsteadyVbc(t,options)
% get FOM pressure based on ROM coefficients
Om_inv = options.grid.Om_inv;

options = set_bc_vectors(t,options);
f       = options.discretization.yM;
dp      = pressure_poisson(f,t,options);
Vbc = - Om_inv.*(options.discretization.G*dp);