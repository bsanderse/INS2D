% project = 'LDC';   % project name used in filenames


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% flow properties
    Re      = 1000;                  % Reynolds number
    visc    = 'laminar';            % laminar or turbulent; 
                                    % influences stress tensor
    nu      = 1/Re;
    regularize = 0;                 %0: no regularization; 1: Leray; 2: C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% domain and mesh
    x1      = 0;
    x2      = 1;
    y1      = 0;
    y2      = 1;

    Nx      = 80;                   % number of volumes in the x-direction
    Ny      = 80;                   % number of volumes in the y-direction

    sx      = 1;                  % stretch factor
    sy      = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time and space discretization

    % steady or unsteady solver
    steady  = 1;         % steady(1) or unsteady(0)

    % spatial accuracy: 2nd or 4th order
    order4  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% solver settings

    % for steady problems or unsteady problems with implicit methods:

    relax                  = 0;    % relaxation parameter to make matrix diagonal more dominant
    
    nonlinear_acc          = 1e-14;
    nonlinear_relacc       = 1e-14;
    nonlinear_maxit        = 10;
        
    nonlinear_Newton       = 1;    % 0: do not compute Jacobian, but approximate iteration matrix with I/dt
                                   % 1: approximate Newton; build Jacobian once at beginning of nonlinear iterations
                                   % 2: full Newton; build Jacobian at each
                                   % iteration
    Jacobian_type          = 1;    % 0: Picard linearization
                                   % 1: Newton linearization
                                   
    % for steady problems only:
    nPicard                = 2;    % in case of Jacobian_type=1, first do nPicard Picard steps                                   
    
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
    
    rtp.show         = 0;          % real time plotting 
    rtp.type         = 'velocity'; % velocity, quiver, vorticity or pressure
    rtp.n            = 1;
    rtp.movie        = 0;
    
%     statistics.write = 1;          % write averages and fluctuations each
%     n steps
%     statistics.n     = 1;
    
    restart.load     = 0;          % start from previous simulation
    restart.folder   = 'results/TCF_100_2x1x1_24x12x12_0';   % folder to be loaded
    restart.file     = 25;         % file number to load
    
    restart.write    = 0;          % write restart files 
    restart.n        = 50;         % every restart.n iterations
    
    save_results     = 1;          % create folder with results files and input files
    path_results     = 'results';  % folder where results are stored
    save_file        = 0;          % save all matlab data after program is completed
    save_unsteady    = 0;
    
    cw_output        = 0;          % command window output; 
                                   % 0: output file, 1: local command window;
                                   % 0 only works if save_results = 1
    
    filelen          = 8;          % number of characters for output files
    
    library_path     = '~/Dropbox/work/Programming/libs/'; % own written matlab libraries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%