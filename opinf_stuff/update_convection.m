function C = update_convection(feasible_basis,A_hat,D,A_dot)

N = size(A_dot,1);

A      = A_hat(1:N,:);
A_kron = A_hat(N+1:end,:);

B_T = (A_dot - D*A)';
B_vec = B_T(:);



CONTINUE!