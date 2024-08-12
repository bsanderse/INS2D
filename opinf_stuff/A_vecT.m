function A_vecT = A_vecT(A,T)
% perform multiplication
% sparse(kron(eye(M),A))*T
% matrix-free

[a1,a2] = size(A);
[t1,t2] = size(T);

T_wide = reshape(T,a2,[]);
A_vecT_wide = A*T_wide;

A_vecT = reshape(A_vecT_wide,[],t2);