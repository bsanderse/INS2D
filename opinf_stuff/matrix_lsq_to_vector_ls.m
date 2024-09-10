function [A,b,T] = matrix_lsq_to_vector_ls(D,X_dot)
% translates matrix-valued least-squares problem into linear system Ax=b
% for vector x

% n = 12
% 
% Q = magic(n);
% % D = [1:n^2; ones(n-1,n^2)]
% D = rand(n,2*n^2);
% X_dot = Q*D;

% Q1 = X_dot/D;
% norm(Q-Q1)

% [q1,q2] = size(Q);
q1 = size(X_dot,1);
q2 = size(D,1);
N = q1*q2;
T = zeros(q1,q2,N);
for i = 1:N
    Ti = zeros(N,1);
    Ti(i) = 1;
    Ti = reshape(Ti,q1,q2);
    T(:,:,i) = Ti;
end

P = N;
K = size(D,2);

A = zeros(P);
b = zeros(P,1);
for p=1:P
    for k=1:K
        b(p) = b(p) + X_dot(:,k)'*T(:,:,p)*D(:,k);
    end
    for s=1:P
        for k=1:K
            A(p,s) = A(p,s) + D(:,k)'*T(:,:,s)'*T(:,:,p)*D(:,k);
            % s,p
            % T(:,:,s)'*T(:,:,p)
        end
    end
end

% c = A\b;
% 
% Q2 = zeros(q1,q2);
% for i = 1:N
%     Q2 = Q2 + c(i)*T(:,:,i);
% end
% 
% norm(Q-Q2)
