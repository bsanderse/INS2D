
A = A_hat';
B = A_dot';

X = A_hat'\A_dot'

[U,S,V] = svd(A_hat','econ');

r = 3

Ut = U(:,1:r);
St = S(1:r,1:r)
Vt = V(:,1:r);

At = Ut*St*Vt';
Aplus_t = Vt*(St\Ut');

Xt = Aplus_t*B

Xt2 = At\B

F1A = St*Vt';
F1B = Ut'*A_dot';
Xt3 = F1A\F1B

r
norm(A*Xt - B,"fro")
norm(At*Xt - B, "fro")

hat_N = size(A_hat,1)

norm(At*Aplus_t*At-At)
norm(Aplus_t*At*Aplus_t-Aplus_t)

reduced_operator(X'- Xt')
reduced_operator(Xt'- Xt2')
reduced_operator(Xt'- Xt3')
