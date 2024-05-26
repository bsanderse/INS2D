function A = ambiguity_constraint(N)
% enforce that [H_i]_:,j = [H_j]_:,i for all i =/= j

A = zeros(N^3,N,N,N); 
% could be made more efficient by eliminating 0-rows,
% but that would make the code messier

for i=1:N
    for j=i+1:N
        for k=1:N
            H = zeros(N,N,N);
            H(k,i,j) = 1;
            H(k,j,i) = -1;
            G = H(:,:)';
            A(:,i,j,k) = G(:);
        end
    end
end

A = A(:,:);

A = A';