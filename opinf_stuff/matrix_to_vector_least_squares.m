Q = magic(2)
A = [1:4; ones(1,4)]
A_dot = Q*A

T1 = [1 0; 0 0];
T2 = [0 1; 0 0];
T3 = [0 0; 1 0];
T4 = [0 0; 0 1];

T = zeros(2,2,2);
T(:,:,1) = T1;
T(:,:,2) = T2;
T(:,:,3) = T3;
T(:,:,4) = T4;

X_dot = A_dot;
D = A;
P = 4;
K = size(A,2);

A = zeros(P);
b = zeros(P,1);
for p=1:P
    for k=1:K
        b(p) = b(p) + X_dot(:,k)'*T(:,:,p)*D(:,k);
    end
    for s=1:P
        for k=1:K
            A(p,s) = A(p,s) + D(:,k)'*T(:,:,s)'*T(:,:,p)*D(:,k);
            % s,k
            % T(:,:,s)'*T(:,:,p)
        end
    end
end

c = A\b