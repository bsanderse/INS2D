function [Diff,Conv] = standardOpInf(A_,dt,nu)

A_dot = (A_(:,2:end) - A_(:,1:end-1))/dt;

        A = A_(:,2:end);

        r = size(A,1);
        K = size(A,2);

        A_kron = zeros(r^2,K);
        for j = 1:K
            a_j = A(:,j);
            A_kron(:,j) = kron(a_j,a_j);
        end

        A_hat = [A', A_kron'];

        O = (A_hat\A_dot')';

        Diff = O(:,1:r)/nu;
        Conv = O(:,r+1:end);
