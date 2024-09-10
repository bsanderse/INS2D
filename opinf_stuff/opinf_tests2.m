clear all
% toy problems for op inf
% 1) 2-dim vector problem x_n+1 = x_n + dt*(nu*l*x_n + q*(x_n otimes x_n))

l = 0*[2 0; 0 3];
q = .2*[eye(2) , eye(2)];
nu = 1;
dt = .01;

K = 2000;
x0 = .005*[1;3];

X = zeros(2,K+1);
X(:,1) = x0;

x = x0;
for i =1:K
    x = x + dt*(nu*l*x + q*kron(x,x));
    X(:,1+i) = x;
end

[Diff,Conv] = standardOpInf(X,dt,nu)


