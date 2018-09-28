load('results/meandering_actuator/AD_Re100_N200_dt1e-3.mat','ufine','vfine','pfine');
u_error = uh-ufine;
v_error = vh-vfine;

p_error = p-pfine;

u_error_i(j) = max(abs(u_error))
v_error_i(j) = max(abs(v_error))
u_error_2(j) = sqrt( sum(u_error.^2)/Nu )
v_error_2(j) = sqrt( sum(v_error.^2)/Nv )

p_error_i(j) = max(abs(p_error))
p_error_2(j) = sqrt( sum(p_error.^2)/Np)    
    