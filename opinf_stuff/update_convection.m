function C = update_convection(feasible_basis,A_hat,D,A_dot)

N = size(A_dot,1);

A      = A_hat(1:N,:);
A_kron = A_hat(N+1:end,:);

B_T = (A_dot - D*A)';
B_vec = B_T(:);

% c = A_vecT(A_kron',feasible_basis)\B_vec; fails if B_vec is cvx object
A_vecT_ = A_vecT(A_kron',feasible_basis);

% normal equations approach
% c = (A_vecT_'*A_vecT_)\(A_vecT_'*B_vec);


% SVD approach
% note: to improve efficiency, this SVD could be precomputed/ performed 
% once outside the CVX SDP solve
c = lsq_SVD(A_vecT_,B_vec);








C = vec2conv(feasible_basis*c,N);

