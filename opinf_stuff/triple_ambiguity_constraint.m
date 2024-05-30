function A = triple_ambiguity_constraint(N)
% enforce that for each triple i j k of pairwise distinct indices
% (ordered i < j < k)
% [H_i]jk = [H_k]ji

% A = zeros(N^3,N,N,N);
A = [];


for i=1:N 
    for j=i+1:N 
        for k=j+1:N 
            H = zeros(N,N,N); % (block) row index, block column index, block index
            H(j,k,i) = 1;
            H(j,i,k) = -1;
            
            G = H(:,:)'; % rearrange to convection-like matricized tensor and transpose to allow column-wise vectorization
            % A(:,i,j,k) = G(:);
            % A(:,i,j,k) = H(:);
            A = [A; G(:)'];
        end
    end
end

% A = A(:,:);
% 
% A = A';