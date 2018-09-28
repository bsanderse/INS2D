version 3.1, June 2012

Finite Volume method for Incompressible Navier-Stokes equations

- second order and fourth order (latter not complete for all BC)
- staggered, non-uniform grid
- 2D 
- general boundary conditions, including unsteady. improved kronecker
  products for interpolation.
- body force: actuator surface methods and immersed boundary method
  moving actuators and adapted RK time integration
- turbulence with RANS k-epsilon model
- turbulence with regularization modeling

- steady solver: mixed Picard / Newton; for turbulent flow solution of
  fully coupled system (u,v,p,k,e). this is in general slow and convergence 
  is not guaranteed; better to use unsteady time-stepping to steady state

- unsteady solver: forward Euler, backward Euler (with Picard linearization)
Adams-Bashforth/Crank-Nicolson, one-leg, fully conservative (with or 
without linearization error), explicit and implicit Runge-Kutta
  for turbulent flows only backward Euler available (timestepping to steady state)

- pressure solver: direct or CG. CG is a mex-file linking a Fortran subroutine

results writing:
- a new folder in the results folder is automatically created for each run
- this folder contains a subfolder with a copy of inputfiles/ at the time 
of running, so that all settings are saved
- these files are used to restart a simulation when restart.load=1. in this 
case NO new folder is created in the results folder.


files that should be edited in order to run simulations:
- parameters.m
- boundary_conditions.m
- initialize.m
- force.m
- mesh_generation.m


to be done:
-temperature effects