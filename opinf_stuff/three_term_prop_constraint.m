function A = three_term_prop_constraint(N)
% actually six term property: enforce that for each tupel i j k, the sum
% over H_p is 0 with p all permutations of i j k

%% original
% A = zeros(N^3,N,N,N); 
% % could be made more efficient by eliminating repetitive permutations,
% % but that would make the code messier
% 
% for i=1:N
%     for j=1:N
%         for k=1:N
%             P = perms([i j k]);
%             H = zeros(N,N,N);
%             for m =1:6
%                 p = P(m,:);
%                 H(p(1),p(2),p(3)) = 1; 
%             end
%             A(:,i,j,k) = H(:);
%         end
%     end
% end
% 
% A = A(:,:);
% 
% % A = A'; A is symmetric (actually w.r.t. to any permutation of i j k,
% % obviously)


%% efficient (=sparse)
if N == 1
    N_c = 1;
elseif N == 2
    N_c = 4;
else
    N_c = nchoosek(N,3);
end
A = spalloc(N^3,N_c,N^3);
counter = 1;

for i=1:N
    for j=1:N
        for k=1:N
            P = perms([i j k]);
            % H = zeros(N,N,N);
            H = spalloc(N^3,1,6);
            for m =1:6 % possible improvement: remove repeating permutations
                p = P(m,:);
                % H(p(1),p(2),p(3)) = 1; 
                H = H + one_hot(p(1),p(2),p(3));
            end
            % A(:,i,j,k) = H(:);
            A(:,counter) = H;
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
        vec = sparse(sub2ind([N,N,N],b,c,a),1,1,N^3,1); 
        % we traverse the third-order tensor first in second, then third, 
        % then first dimension

    end

end