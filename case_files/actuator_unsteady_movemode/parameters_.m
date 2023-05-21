% project = 'actuator_unsteady';   % project name used in filenames
% run_multiple = 0;
run_multiple = 1;
% M_list = [10 10];
% M_list = 8*[10 10 10 10];
% M_list = [10 10 10 10];
% M_list = [10 10 10 20];
% base_M = 20;
% base_M = 10;
base_M = 16;
% M_list = base_M;
% M_list = [base_M base_M base_M 2*base_M];
% M_list = [base_M base_M 2*base_M];
% M_list = [base_M base_M base_M];
% M_list = [base_M base_M base_M 2*base_M];
% M_list = kron(M_list, [1 1]);
% 
% M_list = [2 5 10 20 40 80];%ones(1,5);%kron([2 5 10 20 50 100],ones(1,5));
M_list = [2 5 10 20 40 80];%ones(1,5);%kron([2 5 10 20 50 100],ones(1,5));
% M_list = 2;
% M_list = [2 5 10];
% M_list = [80 79 40];
% M_list = 80;
% M_list = 79;
% M_list = 50;
% M_list = 40;
M_list_raw = M_list;
% repetitions = 5;
% M_list = flip(kron(ones(1,repetitions),M_list_raw));

% M_list = [2 5 8];
M_list = kron(M_list, [1 1]);
% M_list = kron([2 5 10 20 40],[1 1]);
% M_list = [40 40];
% M_list = [60 60];
% M_list = 10*ones(6,1);
% M_list = 10*ones(6,1);
% M_list = 10;
% M_list = 20;
% M_list = 60;

% M_list = [2 5 10 20 2 5 10 20];
% M_list = [2 2 5 5 10 10 20 20];
% M_list = kron([2 5 10 20 50 100],ones(1,5));
% M_list = 100*ones(1,5);
% M_list = flip(kron(ones(1,5),[2 5 10 20 50 100]));
% Mbc_list = [2 4 10 20];
% Mbc = Mbc_list(j);
% Mbc = 10;

mesh_list = ones(length(M_list),1);
% changing_snapshotdata = 1;
changing_snapshotdata = 0;
% if mod(j,2)==0 %j>4 %false %j>4
% % if true %mod(j,2)==0 %j>4 %false %j>4
% % % if false %j>4
% % %     suffix = " mc";
% % %     suffix = " CC";
% % %     suffix = " without lifting function";
% %     suffix = " POD";
% % else
% % %     suffix = " mc";
% %     suffix = " Mbc = "+num2str(Mbc);
% % end
% suffix = " proj div ROM"
% suffix = " Mbc = "+num2str(Mbc);
% suffix = "POD";
% dispName = "ROM M ="+M+suffix;
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
    x1      = 0;
    x2      = 10;
    y1      = -2;
    y2      = 2;

    Nx      = 200; %200;                  % number of volumes in the x-direction
    Ny      = 80; %80;                   % number of volumes in the y-direction

    sx      = 1;                  % stretch factor
    sy      = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% force
    % force is to be set in force.m
    Ct = 0.5; % thrust coefficient actuator disk
    D = 1;     % diameter actuator disk
    
    force_unsteady     = 0;%1; % set to 1 if force is time dependent
    
    % immersed boundary method
    ibm     = 0;
    
    % position of body
    x_c     = 2;
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
%     t_end         = 4*pi*4;        % end time
%     t_end         = 4*pi;        % end time
    t_end         = 20;        % end time

    dt            = t_end/800;      % time step (for explicit methods it can be
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
%     method 21 : generic implicit RK, can also be used for ROM 
if mod(j,2) ~= 0
    method        = 20;
    RK = 'RK44';
%     RK = 'FE11';
%     RK            = 'M2S4R4'; %'FE11'; %'M2S4R4'; %'RK44P2';
else
    method = 21;
    RK = 'GL1';
end
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

    rom     = 1;      % set to 1 to use ROM solver
    pro_rom = 0;     % set to 1 if FOM should provide snapshots for ROM
    M      = M_list(j); %20; %50;    % number of modes used
%     Mbc = M;
%     Mbc = 2;
    % the full snapshotdataset can be reduced by taking as index
    % 1:Nskip:Nsnapshots
    t_sample  = t_end;  % part of snapshot matrix used for building SVD
%     t_sample  = 4*pi;  % part of snapshot matrix used for building SVD
    dt_sample = dt; % frequency of snapshots to be used for SVD
    precompute_convection = 1;%mod(j,2);%1-(j>4);% mod(j,2);%0;
    precompute_diffusion  = 1;%mod(j,2);%1-(j>4);% mod(j,2);%0;
    precompute_force      = 1;%mod(j,2);%1-(j>4);% mod(j,2);%0; 
    precompute_obc        = 1;
%     precompute_convection = mod(j,2);%1-(j>4);% mod(j,2);%0;
%     precompute_diffusion  = mod(j,2);%1-(j>4);% mod(j,2);%0;
%     precompute_force      = mod(j,2);%1-(j>4);% mod(j,2);%0;
%     precompute_obc        = mod(j,2);
%     precompute_convection = 0;
%     precompute_diffusion  = 0;
%     precompute_force      = 0;
%     precompute_obc       = 0;

%     snapshot_data = 'results/actuator_unsteady_movemode_1.000e+02_200x80_fomdata/matlab_data.mat';
    snapshot_data = 'results/actuator_unsteady_movemode_1.000e+02_200x80_finetime/matlab_data.mat';
    

    rom_bc = 2; % 0: homogeneous (no-slip, periodic);
                % 1: non-homogeneous, time-independent;
                % 2: non-homogeneous, time-dependent

%     time_discB4pres_elim = j>2;
% 
%                 bc_recons = [0 2 2 3 4 5];
%                 bc_recons = [3 5]; 
%                 bc_recon = bc_recons(j); Mp = M;

%     bc_recon = 5;
    bc_recon = 3;
% %     bc_recon = 5; M=M+1;
%     bc_recons = kron([1 1 1],[3 5]);
%     bc_recons = kron([1 1 1 1 1 1],[3 5]);
%     bc_recon = bc_recons(j);
%     Mps = 4*[5 10 15 20];
%     Mp = Mps(j);
    Mp = M;
%     Mps = [1 2 3 4];
%     Mp = Mps(j);
%     Mbc = Mps(j) + 1;
%     bc_recon = 4; 
%     bc_recon = 2; %3-2*(j>1); % 2-mod(j,2); %(j>4)+1; %2-mod(j,2); %(j>4)+1;
%     bc_recon = 3; %3-2*(j>1); % 2-mod(j,2); %(j>4)+1; %2-mod(j,2); %(j>4)+1;

    Mbc = M;
%      Mbc = 80;
%      Mbc = 2;
%          Mp = Mbc;

%     bc_recon = 2+mod(j,2); %3-2*(j>1); % 2-mod(j,2); %(j>4)+1; %2-mod(j,2); %(j>4)+1; 
                  % 0: unsteady is always computed by solving a poisson eq
                  % -> supposed to be Sanderse time-independent V inhom approach
                  % 1: Vbc is linearly combined of solutions to Mbc predefined righ-hand sides
                  % 2: no lifting function is used
                  % 3: POD-based Vbc approximation
                  % 4: POD ROM with F_ROM_notvelocityonly
                  % 5: new standard ROM (= with projected mass equation)

    if bc_recon == 3
        suffix = " vo Mbc = "+num2str(Mbc);
        name = "vo M= "+num2str(M)+" Mbc= "+num2str(Mbc);
    elseif bc_recon == 5
        suffix = " vp ";
        name = "vp M= "+num2str(M)+" Mbc= "+num2str(Mbc);     
    else
        suffix = " bc recon = " + bc_recon;
    end
    
%     bases_constructions = ["mthesis" "closest" "optimal" "qr"];
%     bases_constructions = ["mthesis" "optimal" "qr"];
%     bases_constructions = ["mthesis" "optimal"];
%     bases_constructions = [ "optimal"  "mthesis"];
%     bases_constructions = ["mthesis"];
% %     bases_constructions = ["qr"];
%     bases_constructions = [bases_constructions bases_constructions ...
%                            bases_constructions bases_constructions ...
%                            bases_constructions bases_constructions];
% %     bases_constructions = [bases_constructions; bases_constructions];
% %     bases_constructions = bases_constructions(:);
%     bases_construction = bases_constructions(j);

    bases_construction = "mthesis";
%     bases_construction = "closest";
%     bases_construction = "optimal";
%     bases_construction = "qr";
 
    % affects: pressure computation in notvelocityonly RHS computation
    % 0: time derivative of mass equation RHS
    % 1: time difference quations of mass equations RHS's
                
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
    
    rtp.show         = 0;          % 1: real time plotting 
    rtp.n            = 10;
    rtp.movie        = 0;          % requires rtp.show = 1
    rtp.moviename    = 'actuator_unsteady_ROM'; % movie name
    rtp.movierate    = 15;         % frame rate (/s); note one frame is taken every rtp.n timesteps
    
%     show_sigmas = j<=1;
    show_sigmas = 0;

    
%     statistics.write = 1;          % write averages and fluctuations each
%     n steps
%     statistics.n     = 1;
    
    restart.load     = 0;          % start from previous simulation
    restart.folder   = 'results/TCF_100_2x1x1_24x12x12_0';   % folder to be loaded
    restart.file     = 25;         % file number to load
    
    restart.write    = 0;          % write restart files 
    restart.n        = 50;         % every restart.n iterations
    
    save_file        = 0;          % save all matlab data after program is completed    
    path_results     = 'results';  % folder where results are stored
    save_results     = 0;          % create folder with results files and input files
    save_unsteady    = 1;          % save unsteady simulation data at each time step (velocity + pressure) - requires save_file=1
    
    cw_output        = 1;          % command window output; 
                                   % 0: output file, 1: local command window;
                                   % 0 only works if save_results = 1
    
    filelen          = 8;          % number of characters for output files
    
    library_path     = '~/Dropbox/work/Programming/libs/'; % own written matlab libraries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% verbosity
energy_verbosity = 0; % compute unrequired stuff
debug_mode = 0; % perform all kinds of consistency checks -> far more expensive!
% equivalence_cheat = 0; % botch enforcing equivalence of velocity-pressure and velocity-only ROM
equivalence_cheat = 0; % botch enforcing equivalence of velocity-pressure and velocity-only ROM
% equivalence_cheats = kron([1 1 1],[0 1]); 
% equivalence_cheat = equivalence_cheats(j);
% equivalence_cheat = 0;

% if bc_recon == 5 && equivalence_cheat == 0
%     M = 2*M;
% end