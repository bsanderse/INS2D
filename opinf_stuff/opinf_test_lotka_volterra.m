clear all
% toy problems for op inf - Lotka-Volterra
% 1) 2-dim vector problem x_n+1 = x_n + dt*(nu*l*x_n + q*(x_n otimes x_n))

alpha = 4;
beta  = 3;
gamma = 2;
delta = 5;

l = 1*[alpha 0; 0 -gamma];
q = 1*[[0 -beta; 0 delta], zeros(2)];
nu = 1;
dt = .01;

K = 1000;
x0 = [1;3];

X = zeros(2,K+1);
X(:,1) = x0;

x = x0;
for i =1:K
    x = x + dt*(nu*l*x + q*kron(x,x));
    X(:,1+i) = x;
end

plot(X(1,:),X(2,:))


[Diff,Conv] = standardOpInf(X,dt,nu)
[Diff,Conv] = standardOpInf2(X,dt,nu)

