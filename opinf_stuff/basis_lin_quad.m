function B = basis_lin_quad(r)
% construct basis of a such that phi(B) is full rank for
% phi(a) = [ a; kron(a,a) ]

B_r = [ zeros(r,1) eye(r)];
B_r(r,:) = B_r(r,:) +1;

r_star = @(r) r+(r+1)*r/2;
B = zeros(r,r_star(r));
% B(:,r_star(r-1)+r:end) = B_r; 

for r_i = 1:r
    r_bar = 1 + r - r_i;
    rows = r_bar:r;
    cols = r_bar:r+1;
    B_i = B_r(rows,cols);
    B(1:r_i,r_star(r_i-1)+1:r_star(r_i)) = B_i;
end