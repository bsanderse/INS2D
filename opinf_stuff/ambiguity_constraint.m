function A = ambiguity_constraint(N)
% enforce that [H_i]_:,j = [H_j]_:,i for all i =/= j

%% original
A = zeros(N^3,N,N,N); 
% could be made more efficient by eliminating 0-rows,
% but that would make the code messier

% for i=1:N
%     for j=i+1:N
%         for k=1:N
%             H = zeros(N,N,N);
%             H(k,i,j) = 1;
%             H(k,j,i) = -1;
%             G = H(:,:)';
%             A(:,i,j,k) = G(:);
%         end
%     end
% end
% 
% A = A(:,:);
% 
% A = A';

%% efficient (=sparse)

if N<3
    A = [];
else
    N_c = N*nchoosek(N,2);
    A = spalloc(N^3,N_c,2*N_c);
    counter = 1;

    for i=1:N
        for j=i+1:N
            for k=1:N
                % H = zeros(N,N,N);
                % H(k,i,j) = 1;
                % H(k,j,i) = -1;
                % G = H(:,:)';
                % A(:,i,j,k) = G(:);
                A(:,counter) = one_hot(k,i,j) - one_hot(k,j,i);
                counter = counter + 1;
            end
        end
    end

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