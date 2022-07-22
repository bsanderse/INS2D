function [basis,trunc_Sigma_sqrd,no_modes] = POD(snapshot_matrix,no_modes)

[U,S,~] = svd(snapshot_matrix,'econ');
         
%% limit number of modes to meaningful ones
% cond_fac = 1e-6;
cond_fac = 1e-10;
Sigma = diag(S);
no_modes2 = sum(abs(Sigma/Sigma(1))>cond_fac);
no_modes = min(no_modes,no_modes2);
%%

basis = U(:,1:no_modes);
trunc_Sigma_sqrd = (S(1:no_modes,1:no_modes))^2;