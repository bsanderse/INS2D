% solve a possibly rank-deficient least squares problem 
% minimize || A*X -B ||_F^2
% via SVD
function X = lsq_SVD(A,B)

% Uhat = [Us' vectorwise_kron(Us)'];

if issparse(A)
    [U,S,V] = svds(A,min(size(A)));
else
    [U,S,V] = svd(A,'econ');
end

% rank = sum(abs(diag(S))>1e-7)
rank = sum(abs(diag(S))>sqrt(eps()))

X = (V(:,1:rank)*(S(1:rank,1:rank)\U(:,1:rank)'))*B;
