function A = block_skewsymm_constraint(N)
% enforce that each submatrix block is skew-symmetric

% A = zeros(N^3,N,N,N);
A = [];
% could be made more efficient by eliminating 0-rows,
% but that would make the code messier


for i=1:N % block row index
    for j=i:N % block column index
        for k=1:N % block index
            H = zeros(N,N,N);
            if i==j % block diagonal entry
                H(i,j,k) = 1;
            else % block non-diagonal entry
                H(i,j,k) = 1;
                H(j,i,k) = 1;
            end
            G = H(:,:)'; % rearrange to convection-like matricized tensor and transpose
            % A(:,i,j,k) = G(:);
            % A(:,i,j,k) = H(:);
            A = [A; G(:)'];
        end
    end
end

% A = A(:,:);
% 
% A = A';