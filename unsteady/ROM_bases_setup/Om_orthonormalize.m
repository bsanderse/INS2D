function [orthonorm_basis,R,P] = Om_orthonormalize(raw_basis,options)

Om = options.grid.Om;
Om_sqrt = sqrt(Om);
Om_sqrt_inv = 1./Om_sqrt;

[pro_orthonorm_basis,R,P] = orthonormalize(Om_sqrt.*raw_basis);
orthonorm_basis = Om_sqrt_inv.*pro_orthonorm_basis;

%testing
% norm(orthonorm_basis*R-raw_basis(:,P))
