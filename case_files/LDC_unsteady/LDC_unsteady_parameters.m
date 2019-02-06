% input file                
% project = 'LDC_unsteady';   % project name used in filenames


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% flow properties
    Re      = 1000;                  % Reynolds number
    visc    = 'laminar';              % laminar or turbulent; 
                                      % influences stress tensor
    nu      = 1/Re;
    regularize = 0; %0: no regularization; 1: Leray; 2: C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% domain and mesh
    global x1 x2 y1 y2;
    x1      = 0;
    x2      = 1;
    y1      = 0;
    y2      = 1;

    Nx      = 100; %mesh_list(j);         % number of volumes in the x-direction
    Ny      = 100;                   % number of volumes in the y-direction

    L_x     = x2-x1;
    L_y     = y2-y1;
    deltax  = L_x/Nx;               % uniform mesh size x-direction                                   
    deltay  = L_y/Ny;               % uniform mesh size y-direction

    sx      = 1;                  % stretch factor
    sy      = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% force
    % force is to be set in casename_force.m
    
    force_unsteady     = 0; % set to 1 if force is time dependent
    
    % immersed boundary method
    ibm     = 0;
    
    % position of body
    x_c     = 0;
    y_c     = 0;
    D       = 1; % diameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reduced order model

    rom    = 1;      % set to 1 to use ROM solver
    M      = 10;     % number of modes used
    % the full snapshotdataset can be reduced by taking as index
    % 1:Nskip:Nsnapshots
    t_sample  = 1;  % part of snapshot matrix used for building SVD
    dt_sample = 0.01; % frequency of snapshots to be used for SVD
    precompute_convection = 0;
    precompute_diffusion  = 0;
    precompute_force      = 0; 

    snapshot_data = 'results/LDC_unsteady_snapshotdata/matlab_data_Re1000_t1.mat';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time and space discretization

    % steady or unsteady solver
    steady  = 0;         % steady(1) or unsteady(0)

    % spatial accuracy: 2nd or 4th order    
    order4  = 0;

    % only for steady problems:

        linearization = 'Newton';  % Newton or Picard linearization
        nPicard       = 6;         % in case of Newton, first do nPicard Picard steps
        accuracy      = 1e-8;
        relax         = 0;

    
    % only for unsteady problems:

        dt            = 0.02;       % time step (for explicit methods it can be
                                   % determined during running with dynamic_dt)
        t_start       = 0;        % start time
        t_end         = 1;         % end time

        CFL           = 1;              
        timestep.set  = 0;         % time step determined in timestep.m, 
                                   % for explicit methods
        timestep.n    = 1;         % determine dt every timestep.n iterations

        % timestepping method

        % method 1 : Forward Euler: explicit convection and diffusion, 1st order
        % method 2 : IMEX: implicit diffusion (Crank-Nicolson),
        %            explicit convection (Adams-Bashforth),
        %            second order for theta=1/2
        % method 3 : Backward Euler: implicit convection and diffusion, 1st order
        % method 4 : Crank-Nicolson: implicit convection and diffusion,
        %            extrapolated Picard for c^(n+1), 2nd order, almost same as IM1
        % method 5 : explicit one leg beta; 2nd order
        % method 61/62 : Implicit Midpoint: fully conservative
        %                (saddlepoint-system, no pressure correction), 
        %                with (61) or without (62) linearization error
        % method 71/72 : Implicit Midpoint (pressure correction),
        %                with (71) or without (72) linearization error
        % method 81/82 : Explicit Runge-Kutta 4, 
        %                with 1 (81) or 4 (82) pressure-correction steps
        % method 9     : Implicit Runge-Kutta 4 (Gauss 4)
        % method 101/102: Explicit Runge-Kutta 2, Heun (102: pc at each stage)
        % method 111/112: Wray's 3rd order Runge-Kutta with (111) or without (112) Crank-Nicolson
        % method 12 : General explicit RK method (specify tableau in time_RK_gen)
        % method 13  : SDIRK 2-stage, 3rd order
        % method 14 : ARK Radau IIA/B
        % method 15 : Gauss 6
        % method 16 : Lob IIIC 3-stage
        % method 17 : DIRK energy-conserving
        % method 18 : Lob IIICE
        % method 19 : Lob IIIA (CN)
        
        method            = 30;
        theta             = 0.5;
        RK                = 'RK44';
        
        method_startup    = 21;
        method_startup_no = 2; % number of velocity fields necessary for start-up
                               % = equal to order of method
        % only method 2: 
            % theta for diffusion:
%             theta   = 0.5;  % theta=0.5 gives Crank-Nicolson
            % coefficients for explicit convection
            % Adams-Bashforth: alfa1=3/2, alfa2=-1/2 
            % Forward Euler alfa1=1, alfa2=0
%             alfa1   = 3/2;
%             alfa2   = -1/2;
        % only method 5:
%             beta    = 0.1; % should be Reynolds dependent
        % only method 61 and 62
            use_Schur = 0; % solve using Schur complement (pressure correction like)
        % Picard (0) or extrapolated Picard (1) (method 61)
%             EP      = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% solver settings

    % pressure
    poisson          = 1; % 1: direct solver, 
                          % 2: CG with ILU (matlab), 
                          % 3: CG mexfile, 
                          % 4: CG with IC, own Matlab impl.
                          % 5: Petsc
    p_initial        = 1; % calculate pressure field compatible
                          % with the velocity field at t=0
    p_add_solve      = 1; % do additional pressure solve to make it same 
                          % order as velocity

    CG_acc           = 1e-8;  % accuracy for CG (if poisson=2,3,4)
    CG_maxit         = 1000;    % maximum number of iterations for CG

    % diffusion (method 2)
    poisson_diffusion = 1; % options like poisson pressure
    
    % accuracy for non-linear solves (method 62, 72, 9)
    nonlinear_acc          = 1e-10;
    nonlinear_relacc       = 1e-14;
    nonlinear_maxit        = 10;
    nonlinear_Newton       = 2;    % 0: do not compute Jacobian, but approximate iteration matrix with I/dt
                                   % 1: approximate Newton; build Jacobian once at beginning of nonlinear iterations
                                   % 2: full Newton; build Jacobian at each
                                   % iteration
    Jacobian_type          = 1;    % 0: Picard linearization
                                   % 1: Newton linearization
         
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
    
    rtp.show         = 1;          % real time plotting 
    rtp.n            = 25;
    rtp.movie        = 0;          % make movie based on the real time plots
    rtp.moviename    = 'LDC_unsteady_ERK_N40'; % movie name
    rtp.movierate    = 15;         % frame rate (/s); note one frame is taken every rtp.n timesteps
        
%     statistics.write = 1;          % write averages and fluctuations each
%     n steps
%     statistics.n     = 1;
    
    restart.load     = 0;          % start from previous simulation
    restart.folder   = 'results/TCF_100_2x1x1_24x12x12_0';   % folder to be loaded
    restart.file     = 25;         % file number to load
    
    restart.write    = 0;          % write restart files 
    restart.n        = 50;         % every restart.n iterations
    
    save_file        = 0;          % save all matlab data after program is completed    
    path_results     = 'results';  % path where results are stored
    save_results     = 0;          % write information during iterations/timesteps
    save_unsteady    = 1;          % save unsteady simulation data at each time step (velocity + pressure) - requires save_file=1
        
    cw_output        = 1;          % command window output; 
                                   % 0: output file, 1: local command window;
    
    filelen          = 8;          % number of characters for output files
    
    library_path     = '~/Dropbox/work/Programming/libs/'; % own written matlab libraries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%