function phi_bc = get_velo_consis_phi_bc(X_bc,weight_matrix,rank_sensitive)

if (nargin<3)
    rank_sensitive = true;
end

pro_phi_bc = X_bc*weight_matrix';

phi_bc = orthonormalize(pro_phi_bc,rank_sensitive);

