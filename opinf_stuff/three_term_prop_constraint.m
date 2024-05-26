function A = three_term_prop_constraint(N)
% actually six term property: enforce that for each tupel i j k, the sum
% over H_p is 0 with p all permutations of i j k

A = zeros(N^3,N,N,N); 
% could be made more efficient by eliminating repetitive permutations,
% but that would make the code messier

for i=1:N
    for j=1:N
        for k=1:N
            P = perms([i j k]);
            H = zeros(N,N,N);
            for m =1:6
                p = P(m,:);
                H(p(1),p(2),p(3)) = 1; 
            end
            A(:,i,j,k) = H(:);
        end
    end
end

A = A(:,:);

% A = A'; A is symmetric (actually w.r.t. to any permutation of i j k,
% obviously)


