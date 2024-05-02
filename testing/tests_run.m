% script to be run for unit testing
% actual tests are described in tests.m

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

%% tests that do not require mock-up case
snapshot_matrix = eye(3);
no_modes = 3; 
options.grid.Om = ones(3,1);
[basis,S] = Omega_POD(snapshot_matrix,no_modes,options);

norm(basis-eye(3))==0 || error("Omega_POD")
norm(S-eye(3)) == 0 || error("Omega_POD")

% vectorwise_kron
A = zeros(3,3,3);
A(1,1,1) =1;
A(2,2,2) = 1;
A(3,3,3) = 1;
A = A(:,:)';
norm(A-vectorwise_kron(eye(3))) == 0 || error("vectorwise_kron)")

%% mock-up case -> tests are specified in tests.m
% main('actuator_unsteady_consistent_test','testing/case_files',1);