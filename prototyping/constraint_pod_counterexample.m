% X = [1 0; 1 1; -1 1];
X = [1 0; 2 1; -1 1]

[U1,S1,V1] = svd(X)
phi1 = U1(:,1)
bad1 = X-phi1*phi1'*X
bae1 = norm(bad1(:,1))+norm(bad1(:,2))

Q1 = [1; 0; 0];
Xt = X-Q1*Q1'*X
[U2,S2,V2] = svd(Xt)
phi2 = U2(:,1)
bad2 = Xt-phi2*phi2'*Xt
bae2 = norm(bad2(:,1))+norm(bad2(:,2))
