function u = IDFT2(u_hat,Nx,Ny,vec_trunc)
% IDFT: return 2D inverse discrete Fourier transform 
% the definition here is such that the transform is orthonormal
% the input u should be of size Nx_hat*Ny_hat

% the 'extension' parameter N is used to map from a low number of modes
% N_hat to a larger number of spatial points N

% the extension vector vec_trunc can be used instead of N to indicate
% which rows are to be selected, this allows one to for example use a
% selection that is not based on frequency but e.g. on power spectral density
% of the signal, or to skip positive or negative frequencies

if (nargin==1)
    u = IDFT( IDFT(u_hat.') .');
elseif (nargin==2)
    u = IDFT(IDFT(u_hat.',Nx).',Nx);
elseif (nargin==3)
    u = IDFT(IDFT(u_hat.',Ny).',Nx);
else
    error('truncation based on vector indices is not finished for 2D');
end
    


end