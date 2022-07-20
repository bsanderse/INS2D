function [basis,trunc_Sigma_sqrd] = POD(snapshot_matrix,no_modes)

[U,S,~] = svd(snapshot_matrix,'econ');
basis = U(:,1:no_modes);
trunc_Sigma_sqrd = (S(1:no_modes,1:no_modes))^2;