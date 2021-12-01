function u = IRDFT2(u_hat_real,Nx,Ny,vec_trunc)
% IRDFT2: return 2D real inverse discrete Fourier transform
% the definition here is such that the transform is orthonormal

% the 'extension' parameter N is used to map from a low number of modes
% N_hat to a larger number of spatial points N

% the extension vector vec_trunc can be used instead of N to indicate
% which rows were selected when computing u_hat_real, this allows one to for example use a
% selection that is not based on frequency but e.g. on power spectral density
% of the signal, or to skip positive or negative frequencies

if (nargin==1)
    u = IRDFT(IRDFT(u_hat_real.') .');
elseif (nargin==2)
    u = IRDFT(IRDFT(u_hat_real.',Nx).',Nx);
elseif (nargin==3)
    u = IRDFT(IRDFT(u_hat_real.',Ny).',Nx);
else
    error('truncation based on vector indices is not finished for 2D');
end
    


end