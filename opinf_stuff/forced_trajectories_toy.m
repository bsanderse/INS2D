clear all 

p = @(x) [x; kron(x,x)];


%% linear stuff
% a0 = [0;1];

% M = magic(2);
% % M = M+M';
% 
% % f = @(a) M*a;
% % f = @(a) a+1;
% f = @(a) a+[1;3];
% 
% A = [];
% a = a0;
% for i = 1:7
%     A = [A p(a)];
%     a = f(a);
% end
% 
% rank(A)
% A

%% time-dependent force

r = 12;
a0 = ones(r,1);

B0 = eye(r);
B = B0;
for s = 1:r
    offset = zeros(r,1);
    offset(s) = 1;
    B = [B, B0 + 1];
end

% B
r_hat_star = size(B,2);
primes_ = primes(r_hat_star^2);
primes_ = primes_(1:r_hat_star);

g = @(a,t) a + B*cos(t*primes_)';

A = [];
a = a0;
for i = 1:r_hat_star
    A = [A p(a)];
    a = g(a,i);
end

rank(A)
r+ r+r*(r-1)/2 % linear terms + qudratic terms + half of the mixed terms
% A

