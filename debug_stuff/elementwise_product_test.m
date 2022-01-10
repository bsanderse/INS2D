profile on
s = sparse(40000,1);
s(1,1:100) = 1;
s.*s;
diag(s)*s;
profile off
profile viewer