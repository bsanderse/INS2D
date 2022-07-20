function orthonorm_basis = orthonormalize(raw_basis)

[Q,~] = qr(raw_basis,0);
rank_ = rank(raw_basis); 
orthonorm_basis = Q(:,1:rank_);