function phi_hom = homogeneous_projection(phi_h,options)

M = size(phi_h,2);  

L = options.discretization.A;
Gx   = options.discretization.Gx;
Gy   = options.discretization.Gy;
G = [Gx;Gy];
M_h = options.discretization.M;

Om_inv = options.grid.Om_inv;

pro_phi_hom = phi_h + (-Om_inv).*(G*(L\(M_h*phi_h)));

phi_hom = Om_orthonormalize(pro_phi_hom,options);

% testing
norm(M_h*phi_h)
norm(M_h*phi_hom)