function orthonorm_basis = Om_orthonormalize(raw_basis,options)

Om = options.grid.Om;
Om_sqrt = sqrt(Om);
Om_sqrt_inv = 1./Om_sqrt;

orthonorm_basis = Om_sqrt_inv.*orthonormalize(Om_sqrt*raw_basis);
