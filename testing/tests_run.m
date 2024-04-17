% script to be run for unit testing
% actual tests are described in tests.m

%% tests that do not require mock-up case
snapshot_matrix = eye(3);
no_modes = 3; 
options.grid.Om = ones(3,1);
[basis,S] = Omega_POD(snapshot_matrix,no_modes,options);

norm(basis-eye(3))==0 || error("Omega_POD")
norm(S-eye(3)) == 0 || error("Omega_POD")

%% mock-up case -> tests are specified in tests.m
% main('actuator_unsteady_consistent_test','testing/case_files',1);