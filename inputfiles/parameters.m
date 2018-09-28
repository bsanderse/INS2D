% input file                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% flow properties
    u_inf   = 1;
    delta   = 1;
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

    Nx      = 64; %mesh_list(j);         % number of volumes in the x-direction
    Ny      = 64;                   % number of volumes in the y-direction

    L_x     = x2-x1;
    L_y     = y2-y1;
    deltax  = L_x/Nx;               % uniform mesh size x-direction                                   
    deltay  = L_y/Ny;               % uniform mesh size y-direction

    sx      = 1;                  % stretch factor
    sy      = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% force
    % force is to be set in force.m
    
    force_unsteady     = 0; % set to 1 if force is time dependent
    
    % immersed boundary method
    ibm     = 0;
    
    % position of body
    x_c     = 0;
    y_c     = 0;
    D       = 1; % diameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time discretization

    % steady or unsteady solver
    steady  = 1;         % steady(1) or unsteady(0)

    order4  = 0;
    alfa    = 3^4;       % richardson extrapolation factor
    beta    = 9/8;       % interpolation factor: (9/8 for 4th order (only effective for order4=1), 1 for 2nd order)
    
    % only for steady problems:

        linearization = 'Newton';  % Newton or Picard linearization
        nPicard       = 6;         % in case of Newton, first do nPicard Picard steps
        accuracy      = 1e-8;
        relax         = 0;

    
    % only for unsteady problems:

        dt            = 0.01;       % time step (for explicit methods it can be
                                   % determined during running with dynamic_dt)
        t             = 0;        % start time
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
        
        method            = 61;
        
        method_startup    = 61;
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
    poisson          = 2; % 1: direct solver, 
                          % 2: CG with ILU (matlab), 
                          % 3: CG mexfile, 
                          % 4: CG with IC, own Matlab impl.
                          % 5: Petsc
    p_initial        = 0; % calculate pressure field compatible
                          % with the velocity field at t=0
    p_add_solve      = 1; % do additional pressure solve to make it same 
                          % order as velocity

    CG_acc           = 1e-8;  % accuracy for CG (if poisson=2,3,4)
    CG_maxit         = 1000;    % maximum number of iterations for CG

    % diffusion (method 2)
    poisson_diffusion = 1; % options like poisson pressure
    
    % accuracy for non-linear solves (method 62, 72, 9)
    nonlinear_acc          = 1e-14;
    nonlinear_relacc       = 1e-14;
    nonlinear_maxit        = 10;
    nonlinear_build_matrix = 1;    % for small dt one can approximate
                                   % the matrix with I/dt
    nonlinear_Newton       = 1;    % take full Jacobian
    simplified_Newton      = 0;    % constant matrix during nonlinear iteration
    nonlinear_startingvalues = 0;
                                   
    % location of PETSc-matlab mex files                                    
    petsc_mex        ='~/Software/petsc-3.1-p5/bin/matlab/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% output parameters and visualization
    plotgrid         = 0;          % plot gridlines and pressure points
    
    tecplot.write    = 0;          % write to tecplot file
    tecplot.n        = 1;         % write tecplot files every n
    
    rtp.show         = 1;          % real time plotting 
    rtp.type         = 'velocity'; % velocity, quiver, vorticity or pressure
    rtp.n            = 10;
    
%     statistics.write = 1;          % write averages and fluctuations each
%     n steps
%     statistics.n     = 1;
    
    restart.load     = 0;          % start from previous simulation
    restart.folder   = 'results/TCF_100_2x1x1_24x12x12_0';   % folder to be loaded
    restart.file     = 25;         % file number to load
    
    restart.write    = 0;          % write restart files 
    restart.n        = 50;         % every restart.n iterations
    
    save_file        = 0;          % save all matlab data after program is completed
    
    project          = 'LDC';   % project name used in filenames

    path_results     = 'results';  % path where results are stored
    
    cw_output        = 1;          % command window output; 
                                   % 0: output file, 1: local command window;
    
    filelen          = 8;          % number of characters for output files
    
    library_path     = '~/Dropbox/work/Programming/libs/'; % own written matlab libraries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%