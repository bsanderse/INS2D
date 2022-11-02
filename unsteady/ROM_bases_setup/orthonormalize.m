function [orthonorm_basis,R,P] = orthonormalize(raw_basis,rank_sensitive)

if (nargin<2)
    rank_sensitive = true;
end

% [Q,R_,P] = qr(raw_basis,0);
[Q,R_] = qr(raw_basis,0);
P = -666;

if rank_sensitive
    rank_ = rank(raw_basis);
else
    rank_ = size(raw_basis,2);
end
orthonorm_basis = Q(:,1:rank_);
R = R_(1:rank_,:);

%testing
% norm(orthonorm_basis*R-raw_basis(:,P))