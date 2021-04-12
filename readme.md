# Energy-conserving Finite Volume method for the incompressible Navier-Stokes equations
This is the 2D version. [For 3D, click here](https://github.com/bsanderse/INS3D).

## Quick installation and run
- MATLAB is required, but no specific toolboxes. The code has been tested succesfully under Octave in the past (no guarantees).
- Run one of the existing case files by typing `main('casename')`, e.g. `main('LDC')` to run the lid-driven cavity test case
- In general, in order to run simulations:
create a folder in the case_files directory, preferably by copying an existing folder. The files present in this folder should be at least
  - parameters 
  - IC (initial conditions)
  - BCtype (type of boundary conditions)
  - uBC and vBC (boundary condition values)
  - mesh

  and optionally:
  - pp (postprocessing)
  - rtp (real time plot)
  - force (body force)
  - duBCdt and dvBCdt (time-derivative of boundary conditions)

## Spatial discretization
- energy-conserving convection
- second order and fourth order spatial discretization, [paper](https://www.sciencedirect.com/science/article/pii/S0021999113006670)
- staggered cartesian grid (uniform or non-uniform)
- general boundary conditions, including wall, inflow, periodic, outflow, symmetry
- body forces: actuator methods and immersed boundary method
  (also: moving actuators and adapted RK time integration)
- turbulence with RANS-like models, e.g. algebraic eddy viscosity model and k-epsilon model
- turbulence with regularization modeling

## Temporal discretization
- time integration methods: 
   - generic explicit Runge-Kutta, [paper](https://www.sciencedirect.com/science/article/pii/S0021999111006838)
   - generic implicit Runge-Kutta, [paper](https://www.sciencedirect.com/science/article/pii/S002199911200424X)
   - one-leg beta
   - Adams-Bashforth / Crank-Nicolson
   - for RANS only backward Euler available (timestepping to steady state)
- pressure solver: direct (LU), preconditioned CG, or fast Fourier.
- steady solver solves the entire non-linear saddlepoint system with mixed Picard / Newton
  for RANS solution of fully coupled system (u,v,p,k,e) this is in general slow and convergence 
  is not guaranteed; better to use unsteady time-stepping to steady state

## Reduced-order model
- energy-conserving reduced order modeling (ROM) in POD-Galerkin framework, [paper](https://www.sciencedirect.com/science/article/pii/S0021999120305106)
- pressure-free
- precomputation of operators including projection of boundary vectors

## Results and testsuite

### Writing results
There are different ways to get results out of the code.
1. use the syntax `[V,p,options]=main('testcase');` to get the velocity and pressure at the end of the simulation. Useful when the code is called several times from another script, but returns only velocity and pressure (plus code settings given in `options`).
2. use the `save_` options in the parameter file to determine output writing; in this way a new folder in the results folder is created for each run, for example with all time-dependent data and convergence information. Useful for detailed analysis of a testcase with different parameter settings.
3. use the realtimeplotting (rtp) file (executed after each time step). Useful for example to make movies.
4. use the postprocessing file (executed at the end of the simulation). Useful to generate paper-ready plots for a given testcase.

### Running testsuite
To compare with previously obtained benchmark data, use:
- `run_testsuite('case_name')` or
- `run_testsuite({'case_name1','case_name2'})`  