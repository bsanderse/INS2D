function A = block_skewsymm_constraint(N)
% enforce that each submatrix block is skew-symmetric

% A = zeros(N^3,N,N,N);
% A = [];
A = spalloc(N^3,N*(N^2+N)/2,N^3);
counter = 1;

for i=1:N % block row index
    for j=i:N % block column index
        for k=1:N % block index
            H = zeros(N,N,N);
            if i==j % block diagonal entry
                % H(i,j,k) = 1;
                A(:,counter) = one_hot(i,j,k);
                counter = counter + 1;

            else % block non-diagonal entry
                % H(i,j,k) = 1;
                % H(j,i,k) = 1;
                A(:,counter) = one_hot(i,j,k) + one_hot(j,i,k);
                counter = counter + 1;
            end
            % G = H(:,:)'; % rearrange to convection-like matricized tensor and transpose
            % % A(:,i,j,k) = G(:);
            % % A(:,i,j,k) = H(:);
            % A = [A G(:)];
        end
    end
end

% A = A';

    function vec = one_hot(a,b,c)
        % H = zeros(N,N,N); % (block) row index, block column index, block index
        % H(a,b,c) = 1;
        % G = H(:,:)'; % rearrange to convection-like matricized tensor and transpose to allow column-wise vectorization
        % vec = G(:);

        % same thing but sparse:
        vec = sparse(sub2ind([N,N,N],b,a,c),1,1,N^3,1); 
        % we traverse the third-order tensor first in second, then first, 
        % then third dimension

    end
end