function L2 = weightedL2norm(vec,weight)
% input: vec is a matrix of M*N; weight a vector of length M
% output: vector of length M, L2-norm of vector field for each column of
% vec
% for weight=1, the result is the same as the Matlab function vecnorm(vec)

% M: length of vector
% N: number of vectors, e.g. snapshots
[M,N] = size(vec);
L2    = zeros(N,1);
if (length(weight) ~= M)
    error('weight vector should have length corresponding to first dimension of input vector');
end

for i=1:N
    L2(i) = sqrt(vec(:,i)'*(weight.*vec(:,i)));
end

% alternative, vectorized:
% W = spdiags(weight,0,M,M);
% L2 = sqrt(sum(W*(error_V.*error_V)));
