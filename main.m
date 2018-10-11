%   This m-file contains the code for the 2D incompressible Navier-Stokes
%   equations using a Finite Volume Method and a pressure correction
%   method.
%   - horizontal numbering of volumes
%   - 2nd and 4th order spatial (central) discretization convection and diffusion
%   - general boundary conditions; switch fourth order BC by changing
%   addpath('functions/new') to addpath('functions/verstappen'); check
%   operator_convection_diffusion construction of Duy and Dvx

%   see readme.txt

%   Benjamin Sanderse, September 2018


%% close figures and clean variables

clc;
clear all;
close all;
format compact;
format long;
% warning('off');

global uBC vBC;

tic;

%% select case file

% example case names (see case_files directory):
% LDC, BFS, doublejet

folder_cases = 'case_files';
case_name    = 'LDC_unsteady';


%% add paths
addpath('force/');
addpath('functions/new/');
% addpath('functions/verstappen/');
addpath('operators/');
addpath('postprocessing/');
addpath('solvers/');
addpath('time/');
addpath('ibm/');

% path for inputfiles will be determined based on the value of 'restart'
if (~isempty(strfind(path,'inputfiles')))
    rmpath('inputfiles/');
end



%% loop over different meshes


% determine dt for finest mesh
% mesh_list = [40 80];% [160 80 40 20 10];% 20 40 80 160];
% mesh_list = [160 80 40 20 10];
% mesh_list = [5 10 20 40 80];
mesh_list = 64;
% mesh_list = [16 32 64 128];
% dt_list = 0.01; %0.002*pi;
%
for j = 1:length(mesh_list)
    %     Nx = mesh_list(jj);
    
    %     if (jj>1)
    %     dt_list(end+1) = 1/(round(1/(dt_list(end)*2)));
    %     end
    %     dt_list = 0.5*[dt_min dt_min/2 dt_min/10];
    jj = 1;
    % for j = 1:length(dt_list)
    % %     tic;
    %     Nx
    %     dt = dt_list(j)
    %     dt = dt*2
    %     j=1;
    
    % relax = relax_list(j);
    % relax = 0;
    % Re = Re_list(j);
    
    %% load input parameters and constants
    disp(['read input parameters of case ' num2str(case_name)])
    % run('inputfiles/parameters');       % current parameter file
    run([folder_cases '/' case_name '/' case_name '_parameters.m']);
    
    % save into structure 'options'
    accumulate_structure;
    
    % create files and directory for statistics, tecplot, restart, convergence
    % files
    create_files;
    
    if (restart.load == 0)
        
        %     addpath('inputfiles');
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
    % [Fx,Fy] = force(t,options);
    
    
    %% input checking
    input_check;
    
    
    %% start the solver
    fprintf(fcw,['setting-up simulation took ' num2str(toc) '\n']);
    
    fprintf(fcw,'start solver...\n');
    tic
    
    if (steady==1)
        switch visc
            case 'turbulent'
                disp('Steady flow with k-epsilon model, 2nd order');
                solver_steady_ke;
            case 'laminar'
                if (order4==0)
                    if (ibm==0)
                        disp('Steady flow with laminar viscosity model, 2nd order');
                        solver_steady;
                    elseif (ibm==1)
                        disp('Steady flow with laminar viscosity model and immersed boundary method, 2nd order');
                        solver_steady_ibm;
                    else
                        error('wrong value for ibm parameter');
                    end
                elseif (order4==1)
                    disp('Steady flow with laminar viscosity model, 4th order');
%                     solver_steady_4thorder;
                    solver_steady;
                else
                    error('wrong value for order4 parameter');
                end
            otherwise
                error('wrong value for visc parameter');
        end
        
    else
        
        if (strcmp(visc,'turbulent'))
            disp('Unsteady flow with k-eps model, 2nd order');
            solver_unsteady_ke;
        elseif (strcmp(visc,'laminar') || strcmp(visc,'LES'))
            disp('Unsteady flow with laminar or LES model');
            solver_unsteady;
        else
            error('wrong value for visc parameter');
        end
        
        fprintf(fcw,['simulated time: ' num2str(t) '\n']);
        
    end
    
    cpu(j,jj) = toc;
    fprintf(fcw,['total elapsed CPU time: ' num2str(toc) '\n']);
    
    if (poisson==5)
        close(PS);
    end
    
    %% post-processing
    fprintf(fcw,'post-processing...\n');
    post_processing;
    
    % save all data to a matlab file
    if (save_file == 1)
        save(file_mat);
    end
    
    clear n Fx Fy
    % clear nonlinear_its
    %
end
