function [Diff,Conv] = standardOpInf2(A_,dt,nu)

A_dot = (A_(:,2:end) - A_(:,1:end-1))/dt;

        A = A_(:,1:end-1);
%         A = A_(:,2:end);

        r = size(A,1);
        K = size(A,2);

        [~,~,u] = reduced_coordinates(r);

        A_kron = zeros(numel(u),K);
        for j = 1:K
            a_j = A(:,j);
            kron_a_j = kron(a_j,a_j);
            A_kron(:,j) = kron_a_j(u);
        end

        A_hat = [A', A_kron'];

        O = (A_hat\A_dot')';

        Diff = O(:,1:r)/nu;
        Conv = O(:,r+1:end);