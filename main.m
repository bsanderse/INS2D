function [V,p,options] = main(case_name,folder_cases,file_format,param_file_ending)

%   This m-file contains the code for the 2D incompressible Navier-Stokes
%   equations using a Finite Volume Method and a pressure correction
%   method.
%   - horizontal numbering of volumes
%   - 2nd and 4th order spatial (central) discretization convection and diffusion
%   - general boundary conditions; different fourth order BC ('verstappen') can be used by changing
%   addpath('spatial/boundaryconditions/proposed') to addpath('spatial/boundaryconditions/verstappen'); check
%   operator_convection_diffusion construction of Duy and Dvx

%   see readme.txt

%   Benjamin Sanderse, September 2018 - April 2019

if (nargin<1)
    error('please provide an input file');
end

if (nargin==1) ...
        % || sum(folder_cases==0) % botch
    folder_cases = 'case_files'; % default folder name
elseif nargin==2
    file_format = 0;             % default file format
elseif nargin ==3
    param_file_ending = "";      % default parameter file name ending
end

% botch that allows to keep default value for folder_cases while changing
% later inputs
if sum(folder_cases==0) 
       folder_cases = 'case_files'; % default folder name 
end

%% close figures and clean variables

% clc;
% clear vars;
close all;
format compact;
format long;
warning('off','MATLAB:rmpath:DirNotFound');
 

% declare boundary conditions via globals
global uBC vBC dudtBC dvdtBC;

% start timer
tic;


%% add folders to path
addpath('bodyforce/');
addpath('ibm/');
addpath('libs/');
addpath('postprocessing/');
addpath('preprocessing/');
addpath('solvers/pressure/');
addpath('spatial/');
addpath('spatial/operators/');
addpath('spatial/boundaryconditions/');
addpath('spatial/boundaryconditions/proposed/');
addpath('spatial/ROM/');
addpath('spatial/ROM/DEIM');
addpath('steady/');
addpath('unsteady/');
addpath('unsteady/OpInf_ROM_src');
addpath('unsteady/ROM_bases_setup');
addpath('testsuite/');
addpath('results/');

addpath('unsteady/ROM_bases_setup/');
addpath('opinf_stuff/');

% path for inputfiles will be determined based on the value of 'restart'
% if (~isempty(strfind(path,'inputfiles')))
%     rmpath('inputfiles/');
% end


%% load input parameters and constants
disp(['reading input parameters of case: ' num2str(case_name)]);
j = 1; % simulation index counter

if file_format==1
    param_file = [folder_cases '/' case_name '/' 'parameters_' ] +  param_file_ending + ['.m'];
else
    param_file = [folder_cases '/' case_name '/' case_name '_parameters'] + param_file_ending + ['.m'];
end
run(param_file);

% check if multiple simulations should be run
if (~exist('run_multiple','var') || run_multiple == 0)
    Nsim = 1;
else
    Nsim = length(mesh_list);
end
    
% loop over multiple simulations (e.g. different meshes or time steps)
for j=1:Nsim
    
    
    if (j>1)
        % run parameter file again in case we are doing a mesh or parametric study
        run(param_file);
    end
    
    % save into a structure called 'options'
    accumulate_structure;
    
    % create files and directory for statistics, tecplot, restart, convergence
    % files
    create_files;
    
    % remove other cases from the path to prevent them from running
    rmpath(genpath(folder_cases));
    
    if (restart.load == 0)
        
        addpath([folder_cases '/' case_name]);
        
    elseif (restart.load ==1)
        
        fprintf(fcw,['using parameter file from ' restart.folder '\n']);
        addpath([restart.folder '/inputfiles/']);
        % copy restart settings because they will be overwritten
        restart_temp = restart;
        % run parameter file from restart folder
        parameters;
        path_results = restart_temp.folder;
        % copy restart info from current parameter file back
        restart      = restart_temp;
        
    end
    
    
    % add own matlab libraries
    % addpath(library_path);
    
    % add PETSc path (defined in parameters.m)
    if (options.solversettings.poisson == 5)
        addpath(petsc_mex);
    end
    
    
    % turbulence constants
    switch visc
        case 'keps'
            constants_ke;
    end
    
    %% construct mesh
    % disp('construct mesh...')
    fprintf(fcw,'construct mesh...\n');
    mesh_generation;
    
    
    %% boundary conditions
    % disp('boundary conditions...')
    fprintf(fcw,'boundary conditions...\n');
    boundary_conditions;
    
    
    %% construct operators (matrices) which are time-independent
    % disp('construct operators...');
    fprintf(fcw,'construct operators...\n');
    operators;
    
    %% construct matrices mapping ybc to boundary contribution vectors
    options = bc_mappings(options);    
    
    %% initialization of solution vectors
    % disp('initialization of vectors...');
    fprintf(fcw,'initialization of solution vectors...\n');
    if (restart.load == 0)
        initialize;
    end
    
    %% boundary conditions
    options = set_bc_vectors(t,options);
    
    %% initialize basis for lifting function Vbc
    if (options.rom.rom_bc == 2 && options.rom.bc_recon == 1) % could also be used for rom = 0, pro_rom = 1
% ...        || (options.rom.rom_bc == 2 && options.rom.bc_recon == 3) ... %only in debug mode
% ...        || (options.rom.rom_bc == 2 && options.rom.bc_recon == 2)  
        if file_format == 1
            run([folder_cases '/' case_name '/' 'init_unsteady_Vbc_.m']);
        else
            run([folder_cases '/' case_name '/' case_name '_init_unsteady_Vbc.m']);
        end
    end
        
    
    
    %% construct body force or immersed boundary method
    % the body force is called in the residual routines e.g. F.m
    % here we create a handle to bodyforce file, if exists
    if file_format == 1
            file_name_force =  'force_';
    else
            file_name_force = [options.case.project '_force'];
    end
    % file_name_force = [options.case.project '_force'];
    % file_name_force =  'force_';
    if (exist(file_name_force,'file'))
        % create function handle with name bodyforce
        options.force.isforce   = 1;
        options.force.bodyforce = str2func(file_name_force);
        % steady force can be precomputed once:
        if (options.force.force_unsteady==0)
            [options.force.Fx,options.force.Fy] = force(V_start,t,options,0);
        end

    else
        options.force.isforce   = 0;
    end
    
    
    %% input checking
    input_check;
    
    
    %% start the solver
    fprintf(fcw,['setting-up simulation took ' num2str(toc) ' seconds \n']);
    
    fprintf(fcw,'start solver...\n');
    tic
    
    if options.case.test % test mode for unit testing
        disp('RUNNING IN TEST MODE')
    else
        if (options.case.steady==1) % steady

            switch options.case.visc
                case 'keps'
                    fprintf(fcw,'Steady flow with k-epsilon model, 2nd order\n');
                    solver_steady_ke;
                case 'laminar'
                    if (options.discretization.order4==0)
                        if (options.ibm.ibm==0)
                            fprintf(fcw,'Steady flow with laminar viscosity model, 2nd order\n');
                            solver_steady;
                        elseif (options.ibm.ibm==1)
                            fprintf(fcw,'Steady flow with laminar viscosity model and immersed boundary method, 2nd order\n');
                            solver_steady_ibm;
                        else
                            error('wrong value for ibm parameter');
                        end
                    elseif (options.discretization.order4==1)
                        fprintf(fcw,'Steady flow with laminar viscosity model, 4th order\n');
                        %                     solver_steady_4thorder;
                        solver_steady;
                    else
                        error('wrong value for order4 parameter');
                    end
                case 'ML'
                    fprintf(fcw,'Steady flow with mixing length, 2nd order\n');
                    solver_steady;
                otherwise
                    error('wrong value for visc parameter');
            end

        else %unsteady

            switch options.case.visc
                case 'keps'
                    fprintf(fcw,'Unsteady flow with k-eps model, 2nd order\n');
                    solver_unsteady_ke;
                case {'laminar','qr','LES','ML'}
                    if (options.rom.rom==0)
                        fprintf(fcw,'Unsteady flow with laminar or LES model\n');
                        solver_unsteady;
                    elseif (options.rom.rom==1)
                        if options.rom.rom_type == "intrusive"
                            fprintf(fcw,['Unsteady flow with reduced order model with ' num2str(options.rom.M) ' modes \n']);
                            solver_unsteady_ROM;
                        elseif options.rom.rom_type == "OpInf" || options.rom.rom_type == "EC-OpInf"
                            fprintf(fcw,['Unsteady flow with OPERATOR INFERENCE reduced order model with ' num2str(options.rom.M) ' modes \n']);
                            solver_unsteady_ROM_OpInf(snapshot_datas,dt_sample,t_sample,options);
                        else
                            error('wrong value for parameter rom_type')
                        end
                    else
                        error('wrong value for parameter options.rom.rom');
                    end

                otherwise
                    error('wrong value for visc parameter');
            end

            fprintf(fcw,['simulated time: ' num2str(t) '\n']);

        end
    end
    
    cpu(j,1) = toc;
    fprintf(fcw,['total elapsed CPU time: ' num2str(toc) '\n']);
    
    if (options.solversettings.poisson==5)
        close(PS);
    end
    
    if (rtp.movie == 1)
        close(writerObj);
    end

    
    %% post-processing
    fprintf(fcw,'post-processing...\n');
%     post_processing;
    if options.rom.rom_type == "intrusive"
        post_processing2;
    end
    
    % save all data to a matlab file
    if (save_file == 1)
        fprintf(fcw,'saving results to Matlab file...\n');
        save(file_mat);
%         cd(path_results)                             % doesn't work
%         param_file = [case_name '_parameters.m'];
%         save(param_file);
%         cd ../..
    end
    
end