function T = feasible_basis(N)
% construct feasible basis for
% convection operator entries 
% block-skew-symmetry constraint and
% triple-ambiguity constraint

% number of columns = N*(N^2-N) // all block-non-diagonal entries
%                     / 2       // forming skew-symmetry pairs
%                     - N choose 3 // triple ambiguity pairs
no_trip = nchoosek(N,3);

% T = [];
T = spalloc(N^3,N*(N^2-N)/2 - no_trip,N*(N^2-N)/2);

counter = 1;

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
        % T = [T one_hot(i,j,k)-one_hot(j,i,k)];
        T(:,counter) = one_hot(i,j,k)-one_hot(j,i,k);
        counter = counter + 1;

        % case k = j -> no triple ambiguity
        k = j;
        % T = [T one_hot(i,j,k)-one_hot(j,i,k)];
        T(:,counter) = one_hot(i,j,k)-one_hot(j,i,k);
        counter = counter + 1;

        for k=j+1:N % so i < j < k
            H = zeros(N,N,N); % (block) row index, block column index, block index
            
            % combine two of the three skew-symmetry pairs
            % tie breaker: smallest involved row indices
            % T = [T one_hot(i,j,k)-one_hot(j,i,k) ...
            %     + (one_hot(i,k,j)-one_hot(k,i,j))];
            T(:,counter) = one_hot(i,j,k)-one_hot(j,i,k) ...
                        + (one_hot(i,k,j)-one_hot(k,i,j));
            counter = counter + 1;

            % also add third skew-symmetric pair
            % T = [T one_hot(j,k,i)-one_hot(k,j,i)];
            T(:,counter) = one_hot(j,k,i)-one_hot(k,j,i);
            counter = counter + 1;
        end
    end
end

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
