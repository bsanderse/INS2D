% input file                
project = 'LDC';   % project name used in filenames


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% flow properties
    Re      = 100;                  % Reynolds number
    visc    = 'laminar';            % laminar or turbulent; 
                                    % influences stress tensor
    nu      = 1/Re;
    regularize = 0;                 %0: no regularization; 1: Leray; 2: C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% domain and mesh
    global x1 x2 y1 y2;
    x1      = 0;
    x2      = 1;
    y1      = 0;
    y2      = 1;

    Nx      = 40; %mesh_list(j);         % number of volumes in the x-direction
    Ny      = 40;                   % number of volumes in the y-direction

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

    % only for steady problems:
    linearization = 'Newton';  % Newton or Picard linearization
    nPicard       = 6;         % in case of Newton, first do nPicard Picard steps
    accuracy      = 1e-8;
    relax         = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% solver settings
   
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
    
    save_results     = 0;          % create folder with results files and input files
    path_results     = 'results';  % folder where results are stored
    save_file        = 0;          % save all matlab data after program is completed
    
    cw_output        = 1;          % command window output; 
                                   % 0: output file, 1: local command window;
                                   % 0 only works if save_results = 1
    
    filelen          = 8;          % number of characters for output files
    
    library_path     = '~/Dropbox/work/Programming/libs/'; % own written matlab libraries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%