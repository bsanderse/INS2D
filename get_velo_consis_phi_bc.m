function phi_bc = get_velo_consis_phi_bc(X_bc,weight_matrix)

pro_phi_bc = X_bc*weight_matrix';

phi_bc = orthonormalize(pro_phi_bc);

