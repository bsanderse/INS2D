% function options = actuator_unsteady_ROM_init_unsteady_Vbc(t,options)

Om_inv__ = options.grid.Om_inv;
% options = set_bc_vectors(0,options); % already set in main
f__       = options.discretization.yM;
dp__      = pressure_poisson(f__,t,options);
Vbc__     = - Om_inv__.*(options.discretization.G*dp__);

options.rom.Vbc0 = Vbc__;
options.rom.yM0  = options.discretization.yM;
options.rom.abc  = @(t) 1;

