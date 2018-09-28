load(['results/LDC/LDC_ERK4_N' num2str(Nx) '_t1.mat'],'ufine','vfine','pfine');



% temporal error obtained by subtracting spatial error (obtained with
% ufine)
u_error = uh-ufine(:);
v_error = vh-vfine(:);

p_error = p - pfine(:);


u_error_i(j,jj) = max(abs(u_error))
u_error_2(j,jj) = sqrt(sum( u_error.^2)/Nu)
v_error_i(j,jj) = max(abs(v_error))
v_error_2(j,jj) = sqrt(sum( v_error.^2)/Nv)
p_error_i(j,jj) = max(abs(p_error))
p_error_2(j,jj) = sqrt(sum( p_error.^2)/Np)