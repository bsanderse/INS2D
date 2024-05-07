% perform core of OpInf: solve least squares problem and decompose solution
% into linear and quadratic operator
function [L,Q,rank] = OpInf_SVD(Uhats,RHSs)

% Uhat = [Us' vectorwise_kron(Us)'];

[U,S,V] = svd(Uhats,'econ');

rank = sum(abs(diag(S))>1e-7)

O = RHSs*(V(:,1:rank)*(S(1:rank,1:rank)\U(:,1:rank)'));

r = size(RHSs,1);

L = O(:,1:r);
Q = O(:,r+1:end);