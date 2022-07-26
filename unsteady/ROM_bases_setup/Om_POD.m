function [basis,weight_matrix,no_modes] = Om_POD(snapshot_matrix,no_modes,options,cond_fac)
% POD basis with respect to the Omega_h weighted inner product

Om = options.grid.Om;
Om_sqrt = sqrt(Om);
Om_sqrt_inv = 1./Om_sqrt;

[trunc_U,trunc_Sigma_sqrd,no_modes] = POD(Om_sqrt.*snapshot_matrix,no_modes,cond_fac);
basis = Om_sqrt_inv.*trunc_U;

weight_matrix = ((Om_sqrt.*trunc_U)*trunc_Sigma_sqrd)'*snapshot_matrix;

%testing
% basis'*(Om.*basis)

