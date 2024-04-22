clear all
% toy problems for op inf

% 1) scalar problem x_n+1 = x_n + dt*(nu*l*x_n + q*x_n^2)

l = 3;
q = 1;
nu = .1;
dt = .01;

K = 20;
x0 = .005;

X = zeros(1,K+1);
X(1) = x0;

x = x0;
for i =1:K
    x = x + dt*(nu*l*x + q*x^2);
    X(1+i) = x;
end

[Diff,Conv] = standardOpInf(X,dt,nu)
[Diff2,Conv2] = standardOpInf2(X,dt,nu)


