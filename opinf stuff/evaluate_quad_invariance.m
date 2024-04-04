clear all
% define a function that evaluates the invariance over all equivalent
% quadratic operators

N = 4;

N_ = N^2;
p_ = primes(N_);

while length(p_)<N
    N_ = 2*N_;
    p_ = primes(N_);
end

p = p_(1:N);

pp = kron(p,p);

pp0 = unique(pp);

M = repmat(pp0',1,N^2) == pp;

% check
M*pp'./unique(pp)' == sum(M,2)

s_ = (1:numel(pp0))*M;


