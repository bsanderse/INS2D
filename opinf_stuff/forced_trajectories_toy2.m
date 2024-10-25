clear all

rng("default")

p = @(x) [x; kron(x,x)];

r = 8;
r_hat = r + r^2
r_hat_star = r+ r+r*(r-1)/2 % linear terms + qudratic terms + half of the mixed terms


dt = .1;
% M = magic(r); % leads to special cases
M = rand(r,r);

% get symmetric negative semi-definite diffusion dummy
Chol = triu(M);
D = -Chol*Chol';

% get skew-symmetric convection dummy
C = M - M';

M = D - C;


a0 = ones(r,1);

a = a0;
A = zeros(r_hat,r_hat);

a2 = a0;
A2 = zeros(r_hat,r_hat);

for i = 1:r_hat
% for i = 1:r_hat_star
    a = a + dt*M*a;
    A(:,i) = p(a);
    % A(:,i) = p(a)/norm(p(a));

    a2 = a2 + dt*M*a2 + dt*stimulating_function(i*dt,r);
    % a2 = a2 + dt*stimulating_function(i*dt,r);
    A2(:,i) = p(a2);
    % A2(:,i) = p(a2)/norm(p(a2));
end

rank(A)
rank(A2)
    