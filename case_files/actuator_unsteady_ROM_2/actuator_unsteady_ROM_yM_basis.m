% unused

function yMs = actuator_unsteady_ROM_yM_basis(x,y,t,options)
% basis of continuity equation right-hand sides

options = set_bc_vectors(0,options);
yMs = options.discretization.yM;

