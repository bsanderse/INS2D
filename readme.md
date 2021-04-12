# Energy-conserving Finite Volume method for Incompressible Navier-Stokes equations in 2D

## Quick installation and run
- MATLAB is required, no specific toolboxes. Code has been tested succesfully under Octave as well (no guarantees).
- Run one of the existing case files by typing `main('casename')`, e.g. `main('LDC')` to run the lid-driven cavity test case
- In general, in order to run simulations:
create a folder in the case_files directory, preferably by copying an existing folder. The files present in this folder should be at least
  - parameters 
  - IC
  - BCtype
  - uBC
  - vBC
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
   - generic explicit RK, [paper](https://www.sciencedirect.com/science/article/pii/S0021999111006838)
   - generic implicit RK, [paper](https://www.sciencedirect.com/science/article/pii/S002199911200424X)
   - one-leg beta
   - Adams-Bashforth / Crank-Nicolson
   - for RANS only backward Euler available (timestepping to steady state)
- pressure solver: direct (LU), preconditioned CG, or fast Fourier.
- steady solver solves the entire non-linear saddlepoint system with mixed Picard / Newton
  for RANS solution of fully coupled system (u,v,p,k,e) this is in general slow and convergence 
  is not guaranteed; better to use unsteady time-stepping to steady state

## Reduced-order model
- energy-conserving reduced order modeling (ROM) in POD-Galerkin framework, [paper](https://www.sciencedirect.com/science/article/pii/S0021999120305106)
- precomputation of operators

## Results and testsuite

### Writing results
- a new folder in the results folder is automatically created for each run

### Running testsuite
To compare with previously obtained benchmark data, use:
- `run_testsuite('case_name')` or
- `run_testsuite({'case_name1','case_name2'})`  

