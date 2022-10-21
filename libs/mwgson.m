function [Q, R] = mwgson(X,Om)
% Modified WEIGHTED Gram-Schmidt orthonormalization (numerical stable version of Gram-Schmidt algorithm)
% which produces the same result as [Q,R]=qr(X,0)
% Written by Mo Chen (sth4nth@gmail.com).
% further modified by Gustaw Cegielski to account for weighted orthonormalization
Om_mat = spdiags(Om,0,size(Om,1),size(Om,1));
[d,n] = size(X);
m = min(d,n);
R = zeros(m,n);
Q = zeros(d,m);
for i = 1:m
    v = X(:,i);
    for j = 1:i-1
        % R(j,i) = Q(:,j)'*v;
        % v = v-R(j,i)*Q(:,j);
        R(j,i) = Q(:,j)'*Om_mat*v; %projection with weighted norm
        v = v-R(j,i)*Q(:,j);
    end
    % R(i,i) = norm(v);
    R(i,i) = sqrt(v'*Om_mat*v);
    Q(:,i) = v/R(i,i);
end
R(:,m+1:n) = Q'*X(:,m+1:n);