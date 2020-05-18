% project = 'actuator';   % project name used in filenames
run_multiple = 1;
% M_list = [2 2 2 5 5 10 10 20 20 40 40 80 80];
M_list = [2 5 10 20 40 80];
% M_list = 10;
mesh_list = ones(length(M_list),1);
% mesh_list = [1 1 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% flow properties
    Re      = 500;                  % Reynolds number
    visc    = 'laminar';            % laminar or turbulent; 
                                    % influences stress tensor
    nu      = 1/Re;
    regularize = 0;                 %0: no regularization; 1: Leray; 2: C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% domain and mesh
    x1      = -4;
    x2      = 8;
    y1      = -2;
    y2      = 2;

    Nx      = 240;                   % number of volumes in the x-direction
    Ny      = 80;                   % number of volumes in the y-direction

    sx      = 1;                    % stretch factor
    sy      = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% force
    % force is to be set in force.m
    Ct = 0.5; % thrust coefficient actuator disk
    D = 1;     % diameter actuator disk
    
    force_unsteady     = 1; % set to 1 if force is time dependent
    
    % immersed boundary method
    ibm     = 0;
    
    % position of body
    x_c     = 0;
    y_c     = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reduced order model

    rom    = 1;      % set to 1 to use ROM solver
    M      = M_list(j);     % number of velocity modes used
    Mp     = M;     % number of pressure modes used
    % the full snapshotdataset can be reduced by taking as index
    % 1:Nskip:Nsnapshots
    t_sample  = 20;  % part of snapshot matrix used for building SVD
    dt_sample = 10/400; % frequency of snapshots to be used for SVD
    precompute_convection = 1;
    precompute_diffusion  = 1;
    precompute_force      = 1; 
    pressure_recovery     = 1; % compute pressure at each time step
    pressure_precompute   = 1; % precompute PPE operator at ROM level
    pressure_mean         = 0; % subtract mean pressure in constructing ROM
    
    process_iteration_FOM = 1; % execute the process_iteration script each time step (requires FOM evaluation) 
    basis_type            = 1; % 0: choose depending on matrix size, 1: SVD, 2: direct, 3: method of snapshots    
    weighted_norm         = 1;

    rom_bc = 1; % 0: homogeneous (no-slip, periodic); 
                % 1: non-homogeneous, time-independent;
                % 2: non-homogeneous, time-dependent
    snapshot_data = 'results/actuator_ROM_unsteady_force/matlab_data.mat';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time and space discretization

    % steady or unsteady solver
    steady  = 0;         % steady(1) or unsteady(0)

    % spatial accuracy: 2nd or 4th order
    order4  = 0;
    
    % only for unsteady problems:
    dt            = 10/400;      % time step (for explicit methods it can be
                               % determined during running with dynamic_dt)
    t_start       = 0;         % start time
    t_end         = 20; %4*pi;        % end time

    CFL           = 1;              
    timestep.set  = 0;         % time step determined in timestep.m, 
                               % for explicit methods
    timestep.n    = 1;         % determine dt every timestep.n iterations

    % timestepping method

    % method 2 : IMEX AB-CN: implicit diffusion (Crank-Nicolson),
    %            explicit convection (Adams-Bashforth),
    %            second order for theta=1/2
    % method 5 : explicit one leg beta; 2nd order
    % method 20 : generic explicit RK, can also be used for ROM
    % method 21 : generic implicit RK, can also be used for ROM    
    method        = 20;
    RK            = 'RK44P2';

%     method_startup    = 20;
%     method_startup_no = 2; % number of velocity fields necessary for start-up
%     theta = 0.5;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% solver settings

    % pressure (FOM)
    poisson          = 1; % 1: direct solver, 
                          % 2: CG with ILU (matlab), 
                          % 3: CG mexfile, 
                          % 4: CG with IC, own Matlab impl.
                          % 5: Petsc
    p_initial        = 1; % calculate pressure field compatible
                          % with the velocity field at t=0
    p_add_solve      = 1; % do additional pressure solve to make it same 
                          % order as velocity



    % for steady problems or unsteady problems with implicit methods:

    relax                  = 0;    % relaxation parameter to make matrix diagonal more dominant
    
    nonlinear_acc          = 1e-12;
    nonlinear_relacc       = 1e-14;
    nonlinear_maxit        = 10;
        
    nonlinear_Newton       = 1;    % 0: do not compute Jacobian, but approximate iteration matrix with I/dt
                                   % 1: approximate Newton; build Jacobian once at beginning of nonlinear iterations
                                   % 2: full Newton; build Jacobian at each
                                   % iteration
    Jacobian_type          = 1;    % 0: Picard linearization
                                   % 1: Newton linearization
                                   
    % for steady problems only:
    nPicard                = 10;    % in case of Jacobian_type=1, first do nPicard Picard steps                                   
    
    % for unsteady problems only:
    nonlinear_startingvalues = 0;  % extrapolate values from last time step to get accurate initial guess
                                   
    % location of PETSc-matlab mex files                                    
    petsc_mex        ='~/Software/petsc-3.1-p5/bin/matlab/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% output parameters and visualization
    plotgrid         = 0;          % plot gridlines and pressure points
    
    tecplot.write    = 0;          % write to tecplot file
    tecplot.n        = 1;         % write tecplot files every n
    
    rtp.show         = 0;          % 1: real time plotting 
    rtp.n            = 10;
    rtp.movie        = 0;
    rtp.moviename    = 'actuator_ROM'; % movie name
    rtp.movierate    = 15;         % frame rate (/s); note one frame is taken every rtp.n timesteps
        
%     statistics.write = 1;          % write averages and fluctuations each
%     n steps
%     statistics.n     = 1;
    
    restart.load     = 0;          % start from previous simulation
    restart.folder   = 'results/TCF_100_2x1x1_24x12x12_0';   % folder to be loaded
    restart.file     = 25;         % file number to load
    
    restart.write    = 0;          % write restart files 
    restart.n        = 50;         % every restart.n iterations
    
    save_results     = 0;          % create folder with results files (convergence information, pressure solve)
    path_results     = 'results';  % folder where results are stored
    % if save_results=1, the following options can be used:
    save_file        = 0;          % save all matlab data after program is completed in matlab_data.mat
    save_unsteady    = 1;          % save unsteady simulation data at each time step (velocity + pressure) - requires save_file=1
    
    cw_output        = 1;          % command window output; 
                                   % 0: output file, 1: local command window;
                                   % 0 only works if save_results = 1
    
    filelen          = 8;          % number of characters for output files
    
    library_path     = '~/Dropbox/work/Programming/libs/'; % own written matlab libraries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%