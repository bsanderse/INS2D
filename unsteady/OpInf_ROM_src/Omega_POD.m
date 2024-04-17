function [basis,S] = Omega_POD(snapshot_matrix,no_modes,options)
% POD basis with respect to the Omega_h weighted inner product

Om = options.grid.Om;
Om_sqrt = sqrt(Om);
Om_sqrt_inv = 1./Om_sqrt;

[U,S,~] = svd(Om_sqrt.*snapshot_matrix,'econ');
basis = Om_sqrt_inv.*U(:,1:no_modes);