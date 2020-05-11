% test the constrained SVD and compare with Xiao,
% Constrained Proper Orthogonal Decomposition based on QR-factorization for aerodynamical shape optimization

N=40;
p=2;
M=5;

test=zeros(N,p);
test(1:N/2,1)=1;
test(N/2+1:end,2)=1;

% constraint matrix:
G = test/norm(test);

X = rand(N,8);

%% BS approach

% 1) construct (I-ee')*V_svd
Xmod1 = X - G*(G'*X);
% 2) take SVD
[W1,S1,Z1] = svd(Xmod1,'econ');
% 3) add e to get final basis, and truncate
Phi1 = [G W1(:,1:M)];


%% Xiao approach
% constraint: G'*X = 0
% p = rank(G)
[Q,R] = qr(G);
Q1 = Q(:,1:p);
Q2 = Q(:,p+1:end);

% get SVD of Q2'*X
Xmod2 = Q2'*X;
[W2,S2,Z2] = svd(Xmod2,'econ');
Phi2 = [Q1 Q2*W2(:,1:M)];


%% compare
% note that there can be a minus sign difference in the modes
min(abs(Phi1-Phi2),abs(Phi1+Phi2))