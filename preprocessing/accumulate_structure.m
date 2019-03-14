%% fill the structure array options

% in this file the structure options (some sort of global that contains the most
% important settings) is built, and common default values are set for
% the parameters.
% if different parameter values are specified in casename_parameters.m,
% then the ones in that file will be used.

options = [];

%% accumulate options
% object=[];
% 
% voi={
%     
%     };
% 
% accumulate_object;

%% case information
object='case';

voi={
    'project',      case_name;...
    'steady',       [];...    % 0: unsteady; 1: steady
    'visc',         [];...    % 0: laminar; 1: turbulent (k-eps)
    'regularize',   0;...   % convective term regularization; 0: no; 1: Leray; 2: C2
    'force_unsteady', 0;...     % 0: steady forcing or no forcing; 1: unsteady forcing
    'ibm', 0; % 0: no immersed boundary method; 1: immersed boundary method
    };

accumulate_object;

%% physical properties
object='fluid';

voi={
        'Re',      [];
    };
    
accumulate_object;

%% grid parameters
object='grid';

voi={
    'Nx',      [];...
    'Ny',      [];...
    'x',       [];...
    'y',       [];...
    'x1',      [];...
    'x2',      [];...
    'y1',      [];...
    'y2',      [];...
    'deltax',  [];...
    'deltay',  [];...
    'sx',  [];...
    'sy',  [];        
    };

accumulate_object;

%% discretization parameters
object = 'discretization';

voi = {
    'order4',   0;... % 0: 2nd order; 1: 4th order in time 
    'alfa',  81;...   % richardson extrapolation factor = 3^4
    'beta',  9/8;   % interpolation factor
    };

accumulate_object;

%% forcing parameters
object = 'force';

voi = {
    'x_c',  0; ... % x-coordinate of body
    'y_c',  0; ... % y-coordinate of body
    'Ct',   0; ... % thrust coefficient for actuator disk computations
    'D',    1; ... % actuator disk diameter
};

accumulate_object;

%% rom parameters
object = 'rom';

voi = {
    'rom', 0; ... % if 1, use reduced order model
    'M', 10; ... % number of modes for reduced order model
    'precompute_convection', 0; ... % precomputed convection matrices
    'precompute_diffusion', 1; ... % precomputed diffusion matrices
    'precompute_force', 0;... % precomputed forcing term
    't_snapshots', 0; ... % snapshots 
    'dt_snapshots', 0; 
    };

accumulate_object;


%% time marching
object = 'time';

voi = {
    't_start', 0;...
    't_end', [];...
    'dt', []; ...
    'RK', []; ...
    'method', 0;...
    'theta', [];... % theta value for implicit theta method
    };

% options = accumulate_object(object,voi,options);
accumulate_object;

%% solver settings
object = 'solversettings';

voi = {
    'poisson',   1; ... % 1: direct solver; 2: CG with ILU (matlab); 3: CG mexfile; 4: own matlab ICG; 5: PETSC
    'p_initial', 1; ... % 1: calculate compatible IC for the pressure
    'p_add_solve', 1; ... % 1: additional pressure solve to make it same order as velocity
    'CG_acc', 1e-8; ... % accuracy for CG (if poisson=2,3,4)
    'CG_maxit', 1000; ...    % maximum number of iterations for CG
    
    % accuracy for non-linear solves (method 62, 72, 9)
    'nonlinear_acc', 1e-14; ...
    'nonlinear_relacc', 1e-14; ...
    'nonlinear_maxit', 10; ...
    'nonlinear_Newton', 1; ...  % 0: do not compute Jacobian, but approximate iteration matrix with I/dt
                                   % 1: approximate Newton; build Jacobian once at beginning of nonlinear iterations
                                   % 2: full Newton; build Jacobian at each
                                   % iteration
    'Jacobian_type', 0; ...    % 0: Picard linearization, 1: Newton linearization
    'nonlinear_startingvalues', 0; ...
        
    % 
    'nPicard', 5; ... % number of Picard steps before switching to Newton when linearization is Newton
    
    %
    'poisson_diffusion', 1;...
    
    % location of PETSc-matlab mex files                                    
%     petsc_mex        ='~/Software/petsc-3.1-p5/bin/matlab/';
    };

accumulate_object;

%% output files
object = 'output';

voi = {
    'tecplot_write', 0; ...      % write to tecplot file
    'tecplot_n',  1; ...         % write tecplot files every n   
    'save_results', 0; ...
    'path_results', 'results'; ...
    'save_file', 0; ...
    'save_unsteady', 0;
};

accumulate_object;

%% visualization settings
object = 'visualization';

voi={
    'plotgrid', 0; ...       % plot gridlines and pressure points    
    'rtp_show', 1; ...         % real time plotting 
    'rtp_type', 'velocity'; ... % velocity, quiver, vorticity or pressure
    'rtp_n', 10;     
    };

accumulate_object;
