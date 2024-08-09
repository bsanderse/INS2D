function T = feasible_basis(N)
% construct feasible basis for
% convection operator entries 
% block-skew-symmetry constraint and
% triple-ambiguity constraint

T = [];

for i=1:N 
    % case j = i -> block diagonal -> do nothing!!!
    % % % j = i;
    % % % for k= 1:N
    % % %     T = [T one_hot(i,j,k)];
    % % % end

    % case j > 1 -> block non-diagonal
    for j=i+1:N
        % case k = i -> no triple ambiguity
        k = i;
        T = [T one_hot(i,j,k)-one_hot(j,i,k)];

        % case k = j -> no triple ambiguity
        k = j;
        T = [T one_hot(i,j,k)-one_hot(j,i,k)];

        for k=j+1:N % so i < j < k
            H = zeros(N,N,N); % (block) row index, block column index, block index
            
            % combine two of the three skew-symmetry pairs
            % tie breaker: smallest involved row indices
            T = [T one_hot(i,j,k)-one_hot(j,i,k) ...
                + (one_hot(i,k,j)-one_hot(k,i,j))];

            % also add third skew-symmetric pair
            T = [T one_hot(j,k,i)-one_hot(k,j,i)];
        end
    end
end

    function vec = one_hot(a,b,c)
        H = zeros(N,N,N); % (block) row index, block column index, block index
        H(a,b,c) = 1;
        G = H(:,:)'; % rearrange to convection-like matricized tensor and transpose to allow column-wise vectorization
        vec = G(:);
    end
end
