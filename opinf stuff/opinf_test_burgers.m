clear all
% toy problems for op inf - 1D Burgers periodic BC

dx = .1;
N = 20;

D = -2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
D(N,1) = 1;
D(1,N) = 1;
D = D/dx^2;

index = @(x) mod(x-1,N)+1;
kron_ind = @(i,j) i+(j-1)*N;

C = zeros(N,N^2);
for i = 1:N
    C(i,kron_ind(index(i-1),index(i-1))) = 1;
    C(i,kron_ind(index(i+1),index(i+1))) = -1;
    C(i,kron_ind(i,index(i-1))) = 1;
    C(i,kron_ind(i,index(i+1))) = -1;
end

C = C/(3*dx);

C = 0*C;
 
nu = 1;
dt = .01;

K = 100;
x0 = sin(2*pi*(1:N)/N)';

X = zeros(N,K+1);
X(:,1) = x0;

x = x0;
for i =1:K
    x = x + dt*(nu*D*x + C*kron(x,x));
    X(:,1+i) = x;
end



[Diff,Conv] = standardOpInf(X,dt,nu);
[Diff2,Conv2] = standardOpInf2(X,dt,nu);


