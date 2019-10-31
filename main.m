function [V,p] = main(case_name,folder_cases)

%   This m-file contains the code for the 2D incompressible Navier-Stokes
%   equations using a Finite Volume Method and a pressure correction
%   method.
%   - horizontal numbering of volumes
%   - 2nd and 4th order spatial (central) discretization convection and diffusion
%   - general boundary conditions; different fourth order BC ('verstappen') can be used by changing
%   addpath('spatial/boundaryconditions/proposed') to addpath('spatial/boundaryconditions/verstappen'); check
%   operator_convection_diffusion construction of Duy and Dvx

%   see readme.txt

%   Benjamin Sanderse, September 2018 - January 2019

if (nargin<1)
    error('please provide an input file');
end

if (nargin==1) 
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
addpath('postprocessing/');
addpath('preprocessing/');
addpath('solvers/pressure/');
addpath('spatial/');
addpath('spatial/operators/');
addpath('spatial/boundaryconditions/');
addpath('spatial/boundaryconditions/proposed/');
addpath('steady/');
addpath('unsteady/');
addpath('testsuite/');

% path for inputfiles will be determined based on the value of 'restart'
% if (~isempty(strfind(path,'inputfiles')))
%     rmpath('inputfiles/');
% end

   
%% load input parameters and constants
disp(['reading input parameters of case: ' num2str(case_name)]);
j = 1; % simulation index counter
run([folder_cases '/' case_name '/' case_name '_parameters.m']);

% check if multiple simulations should be ran
if (~exist('run_multiple','var') || run_multiple == 0)
    Nsim = 1;
else
    Nsim = length(mesh_list);
end
    
% loop over multiple simulations (e.g. different meshes or time steps)
for j=1:Nsim
    
    
    if (j>1)
        % run parameter file again in case we are doing a mesh study
        run([folder_cases '/' case_name '/' case_name '_parameters.m']);
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
    addpath(library_path);
    
    % add PETSc path (defined in parameters.m)
    if (poisson == 5)
        addpath(petsc_mex);
    end
    
    
    % turbulence constants
    if (strcmp(visc,'turbulent'))
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
    
    
    %% initialization of solution vectors
    % disp('initialization of vectors...');
    fprintf(fcw,'initialization of solution vectors...\n');
    if (restart.load == 0)
        initialize;
    end
    
    %% boundary conditions
    options = set_bc_vectors(t,options);
    
    %% construct body force or immersed boundary method
    % is done in the residual routines e.g. F.m
    
    
    %% input checking
    input_check;
    
    
    %% start the solver
    fprintf(fcw,['setting-up simulation took ' num2str(toc) ' seconds \n']);
    
    fprintf(fcw,'start solver...\n');
    tic
    
    if (steady==1)
        switch visc
            case 'turbulent'
                fprintf(fcw,'Steady flow with k-epsilon model, 2nd order\n');                
                solver_steady_ke;
            case 'laminar'
                if (order4==0)
                    if (ibm==0)
                        fprintf(fcw,'Steady flow with laminar viscosity model, 2nd order\n');
                        solver_steady;
                    elseif (ibm==1)
                        fprintf(fcw,'Steady flow with laminar viscosity model and immersed boundary method, 2nd order\n');
                        solver_steady_ibm;
                    else
                        error('wrong value for ibm parameter');
                    end
                elseif (order4==1)
                    fprintf(fcw,'Steady flow with laminar viscosity model, 4th order\n');
                    %                     solver_steady_4thorder;
                    solver_steady;
                else
                    error('wrong value for order4 parameter');
                end
            otherwise
                error('wrong value for visc parameter');
        end
        
    else
        
        switch visc
            case 'turbulent'
                fprintf(fcw,'Unsteady flow with k-eps model, 2nd order\n');
                solver_unsteady_ke;
            case {'laminar','LES'}
                if (rom==0)
                    fprintf(fcw,'Unsteady flow with laminar or LES model\n');
                    solver_unsteady;
                elseif (rom==1)
                    fprintf(fcw,'Unsteady flow with reduced order model\n');
                    solver_unsteady_ROM;
                else
                    error('wrong value for rom parameter');
                end
                    
            otherwise
                error('wrong value for visc parameter');
        end
        
        fprintf(fcw,['simulated time: ' num2str(t) '\n']);
        
    end
    
    cpu(j,1) = toc;
    fprintf(fcw,['total elapsed CPU time: ' num2str(toc) '\n']);
    
    if (poisson==5)
        close(PS);
    end
    
    if (rtp.movie == 1)
        close(writerObj);
    end

    
    %% post-processing
    fprintf(fcw,'post-processing...\n');
    post_processing;
    
    % save all data to a matlab file
    if (save_file == 1)
        fprintf(fcw,'saving results to Matlab file...\n');
        save(file_mat);
    end
    
    
%     clearvars -except folder_cases case_name j Nsim u_error_i v_error_i ...
%         u_error_2 v_error_2 k_error_i p_error_i p_error_2;
end
