version 4.0, January 2019
Benjamin Sanderse
Finite Volume method for Incompressible Navier-Stokes equations

- second order and fourth order (latter not complete for all BC)
- staggered, non-uniform grid
- 2D 
- general boundary conditions, including unsteady. 
- body force: actuator surface methods and immersed boundary method
  moving actuators and adapted RK time integration
- turbulence with RANS k-epsilon model
- turbulence with regularization modeling
- reduced order modeling (ROM) in POD-Galerkin framework

- steady solver: mixed Picard / Newton; for turbulent flow solution of
  fully coupled system (u,v,p,k,e). this is in general slow and convergence 
  is not guaranteed; better to use unsteady time-stepping to steady state

- unsteady solver: forward Euler, backward Euler (with Picard linearization)
Adams-Bashforth/Crank-Nicolson, one-leg, fully conservative (with or 
without linearization error), explicit and implicit Runge-Kutta
  for RANS only backward Euler available (timestepping to steady state)

- pressure solver: direct or CG. CG is a mex-file linking a Fortran subroutine

results writing:
- a new folder in the results folder is automatically created for each run
- this folder contains a subfolder with a copy of inputfiles/ at the time 
of running, so that all settings are saved
- these files are used to restart a simulation when restart.load=1. in this 
case NO new folder is created in the results folder.


in order to run simulations:
create a folder in the case_files directory, preferably by copying an existing folder
the files present in this folder should be at least
* parameters 
* IC
* BCtype
* uBC
* vBC
* mesh
and optionally
* pp (postprocessing)
* rtp (real time plot)
* force (body force)
* duBCdt (for higher order accurate pressure with time-dependent BC)
* dvBCdt


testsuite:
run_testsuite('case_name') or
run_testsuite({'case_name1','case_name2'})  
to compare with previously obtained benchmark data

to be done:
-temperature effects