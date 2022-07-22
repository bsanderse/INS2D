function [orthonorm_basis,R,P] = orthonormalize(raw_basis)

% [Q,R_,P] = qr(raw_basis,0);
[Q,R_] = qr(raw_basis,0);
P = -666;

rank_ = rank(raw_basis); 
orthonorm_basis = Q(:,1:rank_);
R = R_(1:rank_,:);

%testing
% norm(orthonorm_basis*R-raw_basis(:,P))