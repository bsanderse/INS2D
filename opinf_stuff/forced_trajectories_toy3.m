clear all

rng("default")

p = @(x) [x; -kron(x,x)]; % consistent with convection minus sign policiy opinf implementation

r = 2;
r_hat = r + r^2
r_hat_star = r+ r+r*(r-1)/2 % linear terms + qudratic terms + half of the mixed terms


dt = .1;
% M = magic(r); % leads to special cases
M = rand(r,r);

% get symmetric negative semi-definite diffusion dummy
Chol = triu(M);
D = -Chol*Chol';

% get skew-symmetric convection dummy
C = rand(r,r,r);
C = C - permute(C,[2 1 3]);
C = C(:,:);

a0 = ones(r,1);
nt = r_hat;

a = a0;
A = zeros(r_hat,nt+1);

a2 = a0;
A2 = zeros(r_hat,nt+1);

sforces = zeros(r,nt);

A(:,1) = p(a0);
A2(:,1) = p(a0);

for i = 1:r_hat
% for i = 1:r_hat_star
    a = a + dt*[D C]*p(a);
    A(:,i+1) = p(a);
    % A(:,i) = p(a)/norm(p(a));

    sforce = stimulating_function((i-1)*dt,r);
    sforces(:,i) = sforce;
    a2 = a2 + dt*[D C]*p(a2) + dt*sforce;
    % a2 = a2 + dt*stimulating_function(i*dt,r);
    A2(:,i+1) = p(a2);
    % A2(:,i) = p(a2)/norm(p(a2));
end



rank(A)
rank(A2)

%% infer D and C
dot = @(A) (A(:,2:end)-A(:,1:end-1))/dt;

A2_dot = dot(A2(1:r,:));
A2_dot_clean = A2_dot - sforces;

M_opinf = A2_dot_clean/A2(:,1:end-1);
norm([D C]-M_opinf)

% sanity check
A_dot = dot(A(1:r,:));
M_opinf1 = A_dot/A(:,1:end-1);
norm([D C]-M_opinf1)


%% infer M

% A2_linear = A2(1:r,:);
% rank(A2_linear)
% 
% dot = @(A) (A(:,2:end)-A(:,1:end-1))/dt;
% 
% A2_linear_dot = dot(A2_linear);
% A2_dot_linear_clean = A2_linear_dot - sforces;
% 
% M_opinf = A2_dot_linear_clean/A2_linear(:,1:end-1);
% norm(M-M_opinf)
% 
% % sanity check
% A_linear = A(1:r,:);
% A_dot_linear = dot(A_linear);
% M_opinf1 = A_dot_linear/A_linear(:,1:end-1);
% norm(M-M_opinf1)


    