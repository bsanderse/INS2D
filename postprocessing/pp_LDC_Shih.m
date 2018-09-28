u_exact = 8*(xu.^4 - 2*xu.^3 + xu.^2).*(4*yu.^3 - 2*yu);
v_exact = -8*(4*xv.^3-6*xv.^2+2*xv).*(yv.^4-yv.^2);


u_error = uh-u_exact(:);
v_error = vh-v_exact(:);

% p_error = p - pfine(:);


u_error_i(j,jj) = max(abs(u_error))
u_error_2(j,jj) = sqrt(sum( u_error.^2)/Nu)
v_error_i(j,jj) = max(abs(v_error))
v_error_2(j,jj) = sqrt(sum( v_error.^2)/Nv)
% p_error_i(j,jj) = max(abs(p_error))
% p_error_2(j,jj) = sqrt(sum( p_error.^2)/Np)