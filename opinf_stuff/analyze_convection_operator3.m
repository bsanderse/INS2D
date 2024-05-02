[U,S,V] = svd([A;vectorwise_kron(A)]);
semilogy(diag(S)/S(1,1))