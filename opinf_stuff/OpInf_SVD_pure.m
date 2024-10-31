% perform core of OpInf: solve least squares problem and decompose solution
% into linear and quadratic operator
function [O] = OpInf_SVD_pure(Uhats,RHSs)
% Uhats: r x (r+r^2)
% RHSs: r x K


% Uhat = [Us' vectorwise_kron(Us)'];

if issparse(Uhats)
    [U,S,V] = svds(Uhats,min(size(Uhats)));
else
    [U,S,V] = svd(Uhats,'econ');
end

% rank = sum(abs(diag(S))>1e-7)
rank = sum(abs(diag(S))>sqrt(eps()))

O = RHSs*(V(:,1:rank)*(S(1:rank,1:rank)\U(:,1:rank)'));
% O = (V(:,1:rank)*(S(1:rank,1:rank)\U(:,1:rank)'))*RHSs;
