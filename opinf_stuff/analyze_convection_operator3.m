
A = options.rom.A;
[U,S,V] = svd([A;vectorwise_kron(A)],'econ');
semilogy(diag(S)/S(1,1))