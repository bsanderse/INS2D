function A_kron = vectorwise_kron(A)

r = size(A,1);
K = size(A,2);

A_kron = zeros(r^2,K);
for j = 1:K
    a_j = A(:,j);
    A_kron(:,j) = kron(a_j,a_j);
end