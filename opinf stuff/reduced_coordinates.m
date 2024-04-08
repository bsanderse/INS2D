function s_ = reduced_coordinates(N)

A = reshape(1:N^2,N,N);

B = min(A,A');

b = B(:);

u = unique(b)';


M = repmat(u',1,N^2) == b';

s_ = (1:numel(u))*M;