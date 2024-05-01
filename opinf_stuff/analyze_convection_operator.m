M = size(options.rom.B,2);

C_tensor = reshape(Conv_intrusive,M,M,M);


entries_iii = zeros(M,1);
three_term_prop = zeros(M,1);
two_term_prop = zeros(M);
for i=1:M
    entries_iii(i) = C_tensor(i,i,i);
    U = zeros(M,1);
    U(i) = 1;
    three_term_prop(i) = U'*C_tensor(:,:)*kron(U,U);
    for j =1:M
        V = zeros(M,1);
        V(j) = 1;
        two_term_prop(i,j) = U'*C_tensor(:,:)*kron(V,U);
    end
end

entries_iii
three_term_prop
two_term_prop