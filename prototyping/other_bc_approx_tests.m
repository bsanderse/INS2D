for kk = 1:length(t_js)
    yBC = X_bc(:,kk);
    Y_M3(:,kk) = get_yM(options,yBC);
end

X_inhom = Om_inv.*(G*(L\Y_M3));


t = 0; Vbc_t0 = get_unsteadyVbc(t,options);
norm(X_inhom(:,1)-Vbc_t0) % problem already there

f       = options.discretization.yM;
dp      = pressure_poisson(f,t,options);
Vbc77 = - Om_inv.*(options.discretization.G*dp);