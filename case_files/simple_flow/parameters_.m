% project = 'actuator_unsteady';   % project name used in filenames
run_multiple = 0;
M_list = [10 10];
% M_list = [10 10 10 10];
% M_list = [1 2 5 10 20 40];%ones(1,5);%kron([2 5 10 20 50 100],ones(1,5));
% M_list = [2 5 10 20 40];%ones(1,5);%kron([2 5 10 20 50 100],ones(1,5));
% M_list = kron([2 5 10 20 40],[1 1]);
% M_list = [40 40];
% M_list = 1;
% M_list = 20;


% M_list = [2 5 10 20 2 5 10 20];
% M_list = [2 2 5 5 10 10 20 20];
% M_list = kron([2 5 10 20 50 100],ones(1,5));
% M_list = 100*ones(1,5);
% M_list = flip(kron(ones(1,5),[2 5 10 20 50 100]));
% Mbc_list = [2 4 10 20];
% Mbc = Mbc_list(j);
Mbc = 1;

mesh_list = ones(length(M_list),1);
changing_snapshotdata = 0;
if true %false %mod(j,2)==0 %j>4 %false %j>4
% if false %j>4
%     suffix = " mc";
%     suffix = " CC";
%     suffix = " without lifting function";
    suffix = " POD";
else
%     suffix = " mc";
    suffix = " Mbc = "+num2str(Mbc);
end
% suffix = " Mbc = "+num2str(Mbc);
% suffix = "POD";
% dispName = "ROM M ="+M+suffix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% flow properties
    Re      = 100;                  % Reynolds number
%     Re      = 10^100;                  % Reynolds number
    visc    = 'laminar';            % laminar or turbulent; 
                                    % influences stress tensor
    nu      = 1/Re;
    regularize = 0;                 %0: no regularization; 1: Leray; 2: C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% domain and mesh
    x1      = 0;
    x2      = 10;
    y1      = -2;
    y2      = 2;

    Nx      = 20; %200;                  % number of volumes in the x-direction
    Ny      = 8; %80;                   % number of volumes in the y-direction

    sx      = 1;                  % stretch factor
    sy      = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% force
    % force is to be set in force.m
    Ct = 0; % thrust coefficient actuator disk
    D = 0;     % diameter actuator disk
    
    force_unsteady     = 0;%1; % set to 1 if force is time dependent
    
    % immersed boundary method
    ibm     = 0;
    
    % position of body
    x_c     = 0;
    y_c     = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time and space discretization

    % steady or unsteady solver
    steady  = 0;         % steady(1) or unsteady(0)

    % spatial accuracy: 2nd or 4th order
    order4  = 0;

    % only for unsteady problems:
    t_start       = 0;         % start time
%     t_end         = 4*pi*3;        % end time
    t_end         = 12;        % end time
    dt            = t_end/200;      % time step (for explicit methods it can be
%     dt            = t_end/400;      % time step (for explicit methods it can be
%     dt            = 4*pi/100;      % time step (for explicit methods it can be
                               % determined during running with dynamic_dt)


    CFL           = 1;              
    timestep.set  = 0;         % time step determined in timestep.m, 
                               % for explicit methods
    timestep.n    = 1;         % determine dt every timestep.n iterations
    
    % method 2 : IMEX AB-CN: implicit diffusion (Crank-Nicolson),
    %            explicit convection (Adams-Bashforth),
    %            second order for theta=1/2
    % method 5 : explicit one leg beta; 2nd order
    % method 20 : generic explicit RK, can also be used for ROM
    % method 21 : generic implicit RK, can also be used for ROM    
%     method        = 20;
%     RK = 'FE11';
%     RK = 'RK44';
%     RK            = 'M2S4R4'; %'FE11'; %'M2S4R4'; %'RK44P2';
    method = 21;
    RK = 'GL1';
%     RK = 'BE11';
%     RK = 'RIA1';

    % for methods that are not self-starting, e.g. AB-CN or one-leg
    % beta, we need a startup method.
    % a good choice is for example explicit RK
%     method_startup    = 20;
%     method_startup_no = 2; % number of velocity fields necessary for start-up

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
%     theta = 0.5;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reduced order model

    rom     = 0;      % set to 1 to use ROM solver
    pro_rom = 1;     % set to 1 if FOM should provide snapshots for ROM
    M      = M_list(j); %20; %50;    % number of modes used
    % the full snapshotdataset can be reduced by taking as index
    % 1:Nskip:Nsnapshots
    t_sample  = t_end;  % part of snapshot matrix used for building SVD
    dt_sample = dt; % frequency of snapshots to be used for SVD
%     dt_sample = t_end/200; % frequency of snapshots to be used for SVD
    precompute_convection = 1;%mod(j,2);%1-(j>4);% mod(j,2);%0;
    precompute_diffusion  = 1;%mod(j,2);%1-(j>4);% mod(j,2);%0;
    precompute_force      = 1;%mod(j,2);%1-(j>4);% mod(j,2);%0; 
    precompute_obc       = 1;
%     precompute_convection = mod(j,2);%1-(j>4);% mod(j,2);%0;
%     precompute_diffusion  = mod(j,2);%1-(j>4);% mod(j,2);%0;
%     precompute_force      = mod(j,2);%1-(j>4);% mod(j,2);%0;
%     precompute_convection = 0;
%     precompute_diffusion  = 0;
%     precompute_force      = 0;
%     precompute_obc       = 0;

%     snapshot_data = 'results/actuator_unsteady_snapshotdata/matlab_data.mat';
%     snapshot_data = 'results/actuator_unsteady_ROM_1.000e+02_200x80/matlab_data.mat';
%     snapshot_data = 'results/actuator_unsteady_ROM_1.000e+02_20x8_4/matlab_data.mat'; %M2S4R4
%     snapshot_data = 'results/actuator_unsteady_ROM_1.000e+02_20x8_1/matlab_data.mat'; %   FE11
    
%     snapshot_data = 'results/actuator_unsteady_ROM_1.000e+02_20x8_17/matlab_data.mat'; %M2S4R4
%     snapshot_data = 'results/actuator_unsteady_ROM_1.000e+02_200x80_2/matlab_data.mat'; %M2S4R4

%     snapshot_data = 'results/actuator_unsteady_ROM_1.000e+02_200x80_3/matlab_data.mat'; %M2S4R4 2nd pc

%     snapshot_data =
%     'results/actuator_unsteady_ROM_manipulatedBC/matlab_data.mat';
%     %M2S4R4 2nd pc wrong initial condition
%     snapshot_data = 'results/actuator_unsteady_ROM_manipulatedBC_2/matlab_data.mat'; %M2S4R4 2nd pc
%     snapshot_data = 'results/actuator_unsteady_ROM_1.000e+02_20x8_FOMdata/matlab_data.mat';
%     snapshot_data = 'results/actuator_unsteady_ROM_1.000e+02_20x8_gO=1/matlab_data.mat';
%     snapshot_data = 'results/actuator_unsteady_ROM_1.000e+02_20x8_gO=0/matlab_data.mat';
%     snapshot_data = 'results/actuator_unsteady_ROM_1.000e+02_20x8_debuggedJac/matlab_data.mat';
%     snapshot_data = 'results/actuator_unsteady_ROM_1.000e+02_20x8_gO=1/matlab_data.mat';
% snapshot_data = 'results/simple_flow_1.000e+02_20x8/matlab_data.mat';
% snapshot_data = 'results/simple_flow_1.000e+02_20x8_energy_analysis/matlab_data.mat';
% snapshot_data = 'results/simple_flow_1.000e+02_20x8_e_ana/matlab_data.mat';
% snapshot_data = 'results/simple_flow_1.000e+02_20x8_without_obc/matlab_data.mat';
% snapshot_data = 'results/simple_flow_1.000e+02_20x8_1/matlab_data.mat';
% snapshot_data = 'results/simple_flow_1.000e+02_20x8_FOMdata/matlab_data.mat';
snapshot_data = 'results/simple_flow_1.000e+02_20x8_-obc/matlab_data.mat';

    rom_bc = 1; % 0: homogeneous (no-slip, periodic);
                % 1: non-homogeneous, time-independent;
                % 2: non-homogeneous, time-dependent
    bc_recon = 2; %3-2*(j>1); % 2-mod(j,2); %(j>4)+1; %2-mod(j,2); %(j>4)+1;
%     bc_recon = 3; %3-2*(j>1); % 2-mod(j,2); %(j>4)+1; %2-mod(j,2); %(j>4)+1;
%     bc_recon = 2+mod(j,2); %3-2*(j>1); % 2-mod(j,2); %(j>4)+1; %2-mod(j,2); %(j>4)+1; 
                  % 0: unsteady is always computed by solving a poisson eq
                  % 1: Vbc is linearly combined of solutions to Mbc predefined righ-hand sides
                  % 2: no lifting function is used
                
    process_iteration_FOM = 1; % execute the process_iteration script each time step (requires FOM evaluation) 
    basis_type            = 1; % 0: choose depending on matrix size, 1: SVD, 2: direct, 3: method of snapshots    
    weighted_norm         = 1;
    
    pressure_recovery     =  0; % compute pressure at each time step
    pressure_precompute   =  0; % precompute PPE operator at ROM level
    pressure_mean         =  0; % subtract mean pressure in constructing ROM
    
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
    
    
    % for steady problems or unsteady problems with implicit methods:
    relax                  = 0;    % relaxation parameter to make matrix diagonal more dominant
    
%     nonlinear_acc          = 1e-14;
    nonlinear_acc          = 1e-13;
    nonlinear_relacc       = 1e-14;
%     nonlinear_maxit        = 10;
    nonlinear_maxit        = 100;
        
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
    
    rtp.show         = 1;          % 1: real time plotting 
    rtp.n            = 10;
    rtp.movie        = 1;          % requires rtp.show = 1
    rtp.moviename    = 'simple_flow'; % movie name
    rtp.movierate    = 15;         % frame rate (/s); note one frame is taken every rtp.n timesteps
    
    show_sigmas = 1;
    
%     statistics.write = 1;          % write averages and fluctuations each
%     n steps
%     statistics.n     = 1;
    
    restart.load     = 0;          % start from previous simulation
    restart.folder   = 'results/TCF_100_2x1x1_24x12x12_0';   % folder to be loaded
    restart.file     = 25;         % file number to load
    
    restart.write    = 0;          % write restart files 
    restart.n        = 50;         % every restart.n iterations
    
    save_file        = 1;          % save all matlab data after program is completed    
    path_results     = 'results';  % folder where results are stored
    save_results     = 1;          % create folder with results files and input files
    save_unsteady    = 1;          % save unsteady simulation data at each time step (velocity + pressure) - requires save_file=1
    
    cw_output        = 1;          % command window output; 
                                   % 0: output file, 1: local command window;
                                   % 0 only works if save_results = 1
    
    filelen          = 8;          % number of characters for output files
    
    library_path     = '~/Dropbox/work/Programming/libs/'; % own written matlab libraries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verbosity
energy_verbosity = 0; % compute unrequired stuff
debug_mode = 0; % perform all kinds of consistency checks -> far more expensive!