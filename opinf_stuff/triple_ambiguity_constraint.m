function A = triple_ambiguity_constraint(N)
% enforce that for each triple i j k of pairwise distinct indices
% (ordered i < j < k)
% [H_i]jk = [H_k]ji

% A = zeros(N^3,N,N,N);
% A = [];

if N<3
    A = [];
else
    A = spalloc(N^3,nchoosek(N,3),2*nchoosek(N,3));

    counter = 1;

    for i=1:N
        for j=i+1:N
            for k=j+1:N
                H = zeros(N,N,N); % (block) row index, block column index, block index
                % H(j,k,i) = 1;
                % H(j,i,k) = -1;
                %
                % G = H(:,:)'; % rearrange to convection-like matricized tensor and transpose to allow column-wise vectorization
                % % A(:,i,j,k) = G(:);
                % % A(:,i,j,k) = H(:);
                % A = [A; G(:)'];

                A(:,counter) = one_hot(j,k,i) - one_hot(j,i,k);
                counter = counter +1;
            end
        end
    end

    % A = A(:,:);
    %
    % A = A';
end

    function vec = one_hot(a,b,c)
        % H = zeros(N,N,N); % (block) row index, block column index, block index
        % H(a,b,c) = 1;
        % G = H(:,:)'; % rearrange to convection-like matricized tensor and transpose to allow column-wise vectorization
        % vec = G(:);

        % same thing but sparse:
        vec = sparse(sub2ind([N,N,N],b,c,a),1,1,N^3,1); 
        % we traverse the third-order tensor first in second, then third, 
        % then first dimension

    end

end