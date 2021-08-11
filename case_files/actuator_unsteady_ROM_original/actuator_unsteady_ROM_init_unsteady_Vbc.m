% function options = actuator_unsteady_ROM_init_unsteady_Vbc(t,options)

Om_inv__ = options.grid.Om_inv;
% options = set_bc_vectors(0,options); % already set in main
f__       = options.discretization.yM;
dp__      = pressure_poisson(f__,t,options);
Vbc__     = - Om_inv__.*(options.discretization.G*dp__);

Om__ = options.grid.Om;
normali = norm(sqrt(Om__).*Vbc__);

options.rom.Vbc0 = Vbc__/normali;
options.rom.yM0  = options.discretization.yM;
options.rom.abc  = @(t) cos(pi/6*sin(t/2))*normali;

uBC1 = @(x,y,t,options) actuator_unsteady_ROM_uBC1(x,y,t,options);
vBC1 = @(x,y,t,options) actuator_unsteady_ROM_vBC1(x,y,t,options);
dudtBC1 = @(x,y,t,options) 0*x+0*y;
dvdtBC1 = @(x,y,t,options) 0*x+0*y;
options.bc_options1 = set_bc_vectors_basis(t,options,uBC1,vBC1,dudtBC1,dvdtBC1);
options.rom.abc1  = @(t) cos(pi/6*sin(t/2));

uBC2 = @(x,y,t,options) actuator_unsteady_ROM_uBC2(x,y,t,options);
vBC2 = @(x,y,t,options) actuator_unsteady_ROM_vBC2(x,y,t,options);
dudtBC2 = @(x,y,t,options) 0*x+0*y;
dvdtBC2 = @(x,y,t,options) 0*x+0*y;
options.bc_options2 = set_bc_vectors_basis(t,options,uBC2,vBC2,dudtBC2,dvdtBC2);
options.rom.abc2  = @(t) sin(pi/6*sin(t/2));


% yDiffu1 = bc_options1.discretization.yDiffu;
% yDiffv1 = bc_options1.discretization.yDiffv;
% yDiffu2 = bc_options2.discretization.yDiffu;
% yDiffv2 = bc_options2.discretization.yDiffv;
% norm(yDiffu1-yDiffu2)
% norm(yDiffv1-yDiffv2)

% %% BC base vector 1
% actuator_unsteady_ROM_uBC = @(x,y,t,options) actuator_unsteady_ROM_uBC1(x,y,t,options);
% actuator_unsteady_ROM_vBC = @(x,y,t,options) actuator_unsteady_ROM_vBC1(x,y,t,options);
% options = set_bc_vectors(t,options);
% options.rom.yDiffu1 = options.discretization.yDiffu;
% options.rom.yDiffv1 = options.discretization.yDiffv;
% 
% options.rom.abc1  = @(t) cos(pi/6*sin(t/2));
% 
% 
% %% BC base vector 2
% actuator_unsteady_ROM_uBC = @(x,y,t,options) actuator_unsteady_ROM_uBC2(x,y,t,options);
% actuator_unsteady_ROM_vBC = @(x,y,t,options) actuator_unsteady_ROM_vBC2(x,y,t,options);
% options = set_bc_vectors(t,options);
% options.rom.yDiffu2 = options.discretization.yDiffu;
% options.rom.yDiffv2 = options.discretization.yDiffv;
% 
% options.rom.abc2  = @(t) sin(pi/6*sin(t/2));
% 
% %% restore BC
% actuator_unsteady_ROM_uBC = @(x,y,t,options) actuator_unsteady_ROM_uBC(x,y,t,options); %ohne Sinn :)
% actuator_unsteady_ROM_vBC = @(x,y,t,options) actuator_unsteady_ROM_vBC(x,y,t,options);
% options = set_bc_vectors(t,options);