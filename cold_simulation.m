function [V,p,options] = cold_simulation(case_name,folder_cases,file_format)
% similar to main but without actual simulation

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
        || folder_cases==0 % botch
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
addpath('steady/');
addpath('unsteady/');
addpath('testsuite/');
addpath('debug_stuff/');
addpath('debug_stuff/lions');

% path for inputfiles will be determined based on the value of 'restart'
% if (~isempty(strfind(path,'inputfiles')))
%     rmpath('inputfiles/');
% end

if ~exist('file_format')
    file_format = 0;
end

%% load input parameters and constants
disp(['reading input parameters of case: ' num2str(case_name)]);
j = 1; % simulation index counter
if file_format==1
    run([folder_cases '/' case_name '/' 'parameters_.m']);
else
    run([folder_cases '/' case_name '/' case_name '_parameters.m']);    
end

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
        if file_format==1
                    run([folder_cases '/' case_name '/' 'parameters_.m']);
        else
            run([folder_cases '/' case_name '/' case_name '_parameters.m']);
        end
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
%     file_name_force = [options.case.project '_force'];
    file_name_force =  'force_';
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
    
    
%% main part
divide = @(x) x(2:end)./x(1:end-1);

% vis_velo(options.grid.id_V_inbc2,options);
% vis_p(options.grid.id_p_inbc2,options)

[V_delta,matrix,y_inbc] = background_flow2(options);

M_h = options.discretization.M;
y_M = options.discretization.yM;

mass_cons_error(j) = norm(M_h*V_delta+y_M);

vis_velo(V_delta,options);

matrix_norm(j) = svds(matrix,1);
y_inbc_norm(j) = norm(y_inbc);
V_delta_norm(j) = norm(V_delta);

17


%%
% vis_velo(V,options);
% 
% vis_velo(options.grid.id_V_inbc,options);
% [V_delta,M_inbc,y_inbc] = background_flow1(options);
% 
% M_h = options.discretization.M;
% y_M = options.discretization.yM;
% 
% yM_norms(j) = norm(y_M)
% M_inbc_inv_norm3(j) = 1/sqrt(abs(eigs(M_inbc,1,'smallestabs'))) %wrong
% % M_inbc_inv_norm(j) = 1/sqrt(abs(svds(M_inbc*M_inbc',1,'smallest')))
% M_inbc_inv_norm(j) = 1/(svds(M_inbc,1,'smallest'))
% M_inbc_inv_norm2(j) = norm(full(inv(M_inbc)))
% 
% V_delta_norm(j) = norm(M_inbc\y_inbc)
% V_delta_norm2(j) = norm(inv(M_inbc)*y_inbc)
% y_inbc_norm(j) = norm(y_inbc)
% 
% norm(M_h*V+y_M)
% 
% V_hom = V-V_delta;
% 
% V_hom_norm(j) = norm(V_hom)
% V_hom_norm_omega(j) = sqrt(V_hom'*(Om.*V_hom))
% 
% norm_omega(j) = norm(Om)
% 
% norm(M_h*V_hom)

% vis_velo(V_hom,options)
% vis_velo(V_delta,options)
% vis_velo(V,options)
% 
% V_hom2 = V-V_inhom;
% vis_velo(V_hom2,options)
% vis_velo(V_inhom,options)

%%
% delta_norms(j) = norm(V_delta)
% 
% delta_norms_omega(j) = V_delta'*(Om.*V_delta);
% 
% if j == length(facs)
%     figure
%     loglog(facs,delta_norms)
%     loglog(facs,facs.^(1.5)*10)
%     legend('show')
%     
%     delta_norms./facs.^(1.5)
% end

%%
% K_h = options.discretization.K_h;
% I_h = options.discretization.I_h;
% A_h = options.discretization.A_h;
% y_I = options.discretization.y_I;
% y_A = options.discretization.y_A;
% NF = length(y_A);
% 
% ya_norm(j) = norm(y_A)
% AV_delta_norm(j) = norm(A_h*V_delta)
% advec_norm(j) = norm(A_h*V_delta+y_A)
% VKIV_norm(j) = norm(V_hom'*K_h*spdiags(I_h*V_hom,0,NF,NF))
% 
% advec2(j) = norm(A_h*V_inhom+y_A)
% 
% % b_uvu(j) = V_hom'*K_h*diag(I_h*V_hom)*(A_h*V_delta+y_A)
% b_uvu(j) = V_hom'*K_h*spdiags(I_h*V_hom,0,NF,NF)*(A_h*V_delta+y_A)
% 
% b_uvu2(j) = V_hom'*K_h*spdiags(I_h*V_hom,0,NF,NF)*(A_h*V_inhom+y_A)
% 
% % b_vuu(j) = V_hom'*K_h*diag(I_h*V_delta+y_I)*(A_h*V_hom)
% b_vuu(j) = V_hom'*K_h*spdiags(I_h*V_delta+y_I,0,NF,NF)*(A_h*V_hom)
% 
% b_uuu(j) = V_hom'*K_h*spdiags(I_h*V_hom,0,NF,NF)*(A_h*V_hom)
% 
% matrix1 = K_h*spdiags(I_h*V_hom,0,NF,NF)*(A_h);
% max(matrix1+matrix1',[],'all')
% 
% k_(j) = V'*(Om.*V)
% k_hom(j) = V_hom'*(Om.*V_hom)
% % 
% D_h = options.discretization.D_h;
% diff_(j) = V_hom'*(D_h*V_hom)
% 
% if j == length(facs)
%     figure
% 
%     loglog(facs,abs(b_uvu))
%     hold on
%     loglog(facs,abs(b_vuu))
%     loglog(facs,abs(b_uuu))
%     
%     legend('show')
%     17
% end

%% 

% V_inhom = V;
% V_inhom_norm(j) = norm(V_inhom)
% 
% Om = options.grid.Om;

% L = options.discretization.A;
% L(1,:) = 1;
% Gx   = options.discretization.Gx;
% Gy   = options.discretization.Gy;
% G = [Gx;Gy];
% Om = options.grid.Om;
% Om_inv = options.grid.Om_inv;
% NV = options.grid.NV;
% 
% Om_inv_norm2(j) = norm(Om_inv);
% Om_inv_norm(j) = svds(spdiags(Om_inv,0,NV,NV),1);
% G_norm(j) = svds(G,1);
% L_inv_norm(j) = 1/svds(L,1,'smallest');
% L_inv = inv(L);
% L_inv_norm2(j) = svds(L_inv,1);
% PPE_matrix_norm(j) = svds(Om_inv.*(G*L_inv),1);

%%
% inv_D_h = inv(D_h);
% lamda_h(j) = eigs(inv_D_h,1)

% lambda_h(j) = eigs(D_h,1,'smallestabs')

if j == length(facs)
%     figure
% 
%     loglog(facs,abs(lambda_h))
%     hold on
%     
%     legend('show')
    17
end

%%

17


    
end