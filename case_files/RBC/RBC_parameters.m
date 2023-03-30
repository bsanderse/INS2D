% project = 'RBC';   % project name used in filenames
run_multiple = 1;
% run_multiple = 0;
%mesh_list = [32 64 128 256];
mesh_list = [3e4];
%disp(length(mesh_list));
%mesh_list = [5e-2 1e-2 5e-3 1e-3];
%Nsim = length(mesh_list);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% flow properties
    Re      = 0;                  % Reynolds number
    visc    = 'laminar';            % laminar or turbulent; 
                                    % influences stress tensor
    nu      = 1/Re;
    regularize = 0;                 % 0: no regularization; 1: Leray; 2: C2
    
    boussinesq = 'temp';                 % 'none': mass+ momentum; 'temp' mass+momentum+temperature
    % if boussinesq is used, then the value for Re is not used but
    % calculated from Pr and Ra
    Pr = 0.71;                  % Prandtl number
%     Ra = 1e6;                   % Rayleigh number
    Ra = mesh_list(j);	
    Ge = 1;                   % Gebhart number
    incl_dissipation = 0;       % use dissipation term in temperature equation (1=yes,0=no)
    nondim_type = 1;            % see thermal_constants.m: 1 => uref= sqrt(beta*g*Delta T*H), 2=> uref = kappa/H, 3=> uref = sqrt(c*DeltaT)
    fprintf('Ra =');
    disp(Ra);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% domain and mesh
    x1      = 0;
    x2      = 1;
    y1      = 0;
    y2      = 1;

%    Nx=mesh_list(j);
%    Ny=Nx;
    Nx      = 4;                  % number of volumes in the x-direction
    Ny      = 4;                   % number of volumes in the y-direction
%	if(j==2)
%	Nx=256;
%	Ny=Nx;
%	end
    sx      = 1;                  % stretch factor
    sy      = 1;
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

        dt            = 5e-2;       % time step (for explicit methods it can be
%         dt            = 5e-1;       % time step (for explicit methods it can be
                                   % determined during running with dynamic_dt)
%        dt=mesh_list(j); 
        t_start       = 0;        % start time
        t_end         = 30;         % end time

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
        % method 22 : implicit midpoint (energy conserving) for Boussinesq
        %             system
%         method            = 2;
        method            = 20;
        RK                = 'RK44'; % only used when method = 20 or 21
        
        % for methods that are not self-starting, e.g. AB-CN or one-leg
        % beta, we need a startup method.
        % a good choice is for example explicit RK        
        method_startup    = 20;
        method_startup_no = 1; % number of velocity fields necessary for start-up
                               % = equal to order of method
        % parameters for time integration methods:
        % Adams Bashforth - Crank Nicolson (method 2):
            % theta for diffusion:
%             theta   = 0.5;  % theta=0.5 gives Crank-Nicolson
            % coefficients for explicit convection
            % Adams-Bashforth: alfa1=3/2, alfa2=-1/2 
            % Forward Euler alfa1=1, alfa2=0
%             alfa1   = 3/2;
%             alfa2   = -1/2;
        % one-leg beta (method 5):
%             beta    = 0.5; % in fact, this should be Reynolds dependent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% solver settings

    % pressure
    poisson          = 1; % 1: direct solver, 
                          % 2: CG with ILU (matlab), 
                          % 3: CG mexfile, 
                          % 4: CG with IC, own Matlab impl.
                          % 5: Petsc
                          % 6: FFT                                                    
    p_initial        = 1; % calculate pressure field compatible
                          % with the velocity field at t=0
    p_add_solve      = 0; % do additional pressure solve to make it same 
                          % order as velocity

    CG_acc           = 1e-8;  % accuracy for CG (if poisson=2,3,4)
    CG_maxit         = 1000;    % maximum number of iterations for CG
    
    % for steady problems or unsteady problems with implicit methods:

    relax                  = 0;    % relaxation parameter to make matrix diagonal more dominant
    
    nonlinear_acc          = 1e-12;
    nonlinear_relacc       = 1e-14;
    nonlinear_maxit        = 10;
        
    nonlinear_Newton       = 2;    % 0: do not compute Jacobian, but approximate iteration matrix with I/dt
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
    
    rtp.show         = 1;          % real time plotting 
    rtp.type         = 'velocity'; % velocity, quiver, vorticity or pressure
    rtp.n            = 50;
    rtp.movie        = 0;
    rtp.moviename    = 'RBC';
    rtp.movierate    = 10;           % frame rate (/s); note one frame is taken every rtp.n timesteps
    
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
    save_file        = 1;          % save all matlab data after program is completed
    save_unsteady    = 1;
    sampling_start   = 10/dt; 	   % Time to start sampling, (divided by dt) and it is basically an index not time
    sampling_start_time = sampling_start*dt;
    cw_output        = 1;          % command window output; 
                                   % 0: output file, 1: local command window;
                                   % 0 only works if save_results = 1
    
    filelen          = 8;          % number of characters for output files
    
    library_path     = '~/Dropbox/work/Programming/libs/'; % own written matlab libraries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reduced order model

%     rom    = 0;      % set to 1 to use ROM solver
    rom    = 1;      % set to 1 to use ROM solver

    run_multiple = 0;  % set to 1 to avoid loading FOM data multiple times
%     M_list = [2 4 8 16 2 4 8 16];
%     M      = M_list(j);     % number of modes used
    M      = 28;     % number of modes used (Total velocity modes = Nu+Nv; and Nu=Nx*Ny but Nv= (Ny-1)*Nx
                       % therefore M=2*Nx*Ny-Nx
    Mp     = 16;     % number of pressure modes used (only needed if pressure_recovery=1)
    MT     = 16;     % number of temperature modes (MT=MP=Nx*Ny) used (only needed if boussinesq='temp')

    t_sample  = t_end;%t_end; %t_end;  % part of snapshot matrix used for building SVD
%     t_sample  = 10;%t_end; %t_end;  % part of snapshot matrix used for building SVD
    dt_sample = dt; % frequency of snapshots to be used for SVD
%     dt_sample = 100*dt; % frequency of snapshots to be used for SVD
    take_all_snapshots  =   1; % 1 refers to collection of all the samples; t_sample and dt_sample are not used
                               % 0 refers to only a section of it as specified by t_sample and dt_sample
    initialize_with_snapshots = 1; %
    weighted_norm         = 1;
    weighted_norm_T       = 1; % weighted norm for temeprature basis (only needed if boussinesq='temp')
    basis_type            = 1; % 0: choose depending on matrix size, 1: SVD, 2: direct, 3: method of snapshots
    mom_cons              = 0; %j>4;

    rom_bc = 0; % 0: homogeneous (no-slip, periodic); 
                % 1: non-homogeneous, time-independent;
                % 2: non-homogeneous, time-dependent 

    precompute_convection = 1;
    precompute_diffusion  = 1;
    precompute_force      = 0;
    precompute_buoyancy_force      = 1;

    precompute_convectionT = 0;
    precompute_diffusionT  = 0;

    pressure_recovery     = 0;
    pressure_precompute   = 0;

    process_iteration_FOM = 1; % execute the process_iteration script each time step (requires FOM evaluation)

    snapshot_data = 'RBC_data/matlab_data.mat';
%     snapshot_data = '/export/scratch1/krishan/INS2D_scratch/FOMdata/matlab_data3e4_4x4.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
