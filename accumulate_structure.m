%% fill the structure array options

% in this file the structure options (some sort of global that contains the most
% important settings) is built, and common default values are set for
% the parameters.
% if different parameter values are specified in casename_parameters.m,
% then the ones in that file will be used.

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
    'project',      [];...
    'steady',       [];...
    'visc',         [];    % 0: laminar; 1: turbulent (k-eps)
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
    'order4',   0; % 0: 2nd order; 1: 4th order in time 
    };

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
    'nonlinear_build_matrix', 1; ...    % for small dt one can approximate
                                   % the matrix with I/dt
    'nonlinear_Newton', 1; ...    % take full Jacobian
    'simplified_Newton', 0; ...    % constant matrix during nonlinear iteration
    'nonlinear_startingvalues', 0;
                                   
    % location of PETSc-matlab mex files                                    
%     petsc_mex        ='~/Software/petsc-3.1-p5/bin/matlab/';
    };

accumulate_object;


%% visualization settings
object = 'visualization';

voi={
    'plotgrid', 0; ...       % plot gridlines and pressure points
    
    'tecplot_write', 0; ...          % write to tecplot file
    'tecplot_n',  1; ...         % write tecplot files every n
    
    'rtp_show', 1; ...         % real time plotting 
    'rtp_type', 'velocity'; % velocity, quiver, vorticity or pressure
    'rtp_n', 10; 
    
    };

accumulate_object;